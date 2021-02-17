% colors is (3 vertices) (3 rgb) (nT triangles)
function [extra, energy, colors, gradient] = linearComputeEnergy(img, mesh, n1D, salmap)
    extra = {};
    X = mesh.X; T = mesh.T; nT = size(T,1);
    % generate sample locations in barycentric coords
    [ws, interiorInds] = getBarycentricSamplingWeights(n1D); 
    samplePoints = getSamplePointsFromBarycentricWeights(ws, mesh);
    n = size(ws,1); 
    dA = mesh.triAreas/n; % differential unit of area
    
    % perform sampling
    f_triangle = permute(double(sampleImage(img, samplePoints)),[2 1 3]); % nT, n, 3
    
    % f_triangle: nT  n     3
    % ws:             n  3
    % compute Lj: nT     3  3 
    % triangles, basis funcs, rgb colors
    Lj = squeeze(sum(reshape(f_triangle,nT,n,1,3) .* reshape(ws,1,n,3,1),2)).*dA;
    
    [~, Ki] = loadMassMatrix(1);
    % colors: nT vert rgb
    % Ki     [1   3  3  1 ]
    % Lj_1 = [nT  1  3  3 ]
    % A =    [nT          ]
    extra.colors_fem = squeeze(sum(reshape(Ki,1, 3, 3, 1).*reshape(Lj,nT,1,3,3)./(2*mesh.triAreas), 3));
    
    % Least Squares color per triangle. 
    % In comparison, FEM method doesn't work as well because inverse mass matrix assumes a good numerical integral over each triangle. 
    % Turns out uniform samples for integrating on triangles isn't that accurate.
    [colors, extra.colorsAlt] = linearGetColorsFromSamples(mesh, samplePoints, f_triangle, interiorInds);
    
    % compute actual energy.
    % ws:         [   n  3  ]
    % colors:     [nT    3 3]
    linearColor = squeeze(sum(reshape(ws,1,n,3,1).*reshape(colors,nT,1,3,3),3));
    extra.perTriangleRGBError = dA.*squeeze(vecnorm(linearColor - f_triangle,2,2).^2);
    energy = sum(extra.perTriangleRGBError,'all');
    
    if nargout >= 4
        %% compute gradient
        Lj_1 = permute(Lj,[1 3 2]); % permuted for convenience: nT, rgb, phi
        
        % constantGradients: [nT 6=(2(xy) 3(verts))]
        constantGradients = mesh.dAdt./mesh.triAreas;
        
        % constantGradients: [nT             6 = (2(xy) 3(verts)) ]
        % dphi_dt: triangles, samples, phi, v
        % ws                   [1  n  1         1      3(phi)]  
        % constantGradients_1  [nT 1  3(verts)  2(xy)  1     ]
        % After optimization the code kind of obfuscates the math. verts and phi are supposed to switch though. 
        constantGradients_1 = permute(reshape(constantGradients,nT,1,2,3),[1 2 4 3]); 
        dphi_dt = reshape(squeeze(constantGradients_1.*reshape(ws,1,n,1,1,3)), nT, n, 3, 6);
        
        % vn is (nT, n, 3, 6). 3 for number of edges. 6 for number of velocity values.
        vndl = sampleVdotN_dl(mesh, n1D);

        % generate edge samples of f.
        % edgeSamplePoints: [nT, n, 3(edges), 2(xy)]
        edgeSamplePoints = getTriEdgeSamplePoints(mesh,n1D);
        % f_tri_edges:      [nT, n, 3(edges), 3(rgb)]
        f_tri_edges = double(sampleImage(img, edgeSamplePoints));

        % generate edge samples of phi
        % phi_edges_oneTri: (samples per edge) (3 edges) (3 phis)
        phi_edges_oneTri = getEdgeBarycentricSamplingWeights(n1D);
        
        % phi_edges_oneTri: [      n1D    3(edges)                 3(phis)]
        % vndl:             [nT    n1D    3(edges)   6(v)                 ]
        % f_tri_edges:      [nT    n1D    3(edges)          3(rgb)        ]
        % int_fvn_dl:       [nT)                     6(v)   3(rgb) 3(phis)]
        int_fvnphi_dl = squeeze(sum(reshape(phi_edges_oneTri,1,n1D,3,1,1,3).*vndl.*reshape(f_tri_edges,nT,n1D,3,1,3,1),[2 3]));
        
        % build: dLdt
        % f_triangle_1:  [nT n   3(rgb)         ]
        % dphi_dt_1:     [nT n 6         3(phi) ]
        % int_fvnphi_dl: [nT   6 3(rgb)  3(phi) ]
        % dLdt:          [nT   6 3(rgb)  3(phi) ]
        dphi_dt_1 = reshape(permute(dphi_dt,[1 2 4 3]),nT,n,6,1,3);
        f_triangle_1 = reshape(f_triangle,nT,n,1,3,1);
        % Not entirely sure why this area component needs to have a negative sign. Might have forgotten a sign earlier.
        area_component = -squeeze(sum(dphi_dt_1.*f_triangle_1,2)).*dA; 
        dLdt = area_component + int_fvnphi_dl;
        
        % build: dcdt = Ki/2 * (dLdt/mesh.triAreas - Lj.*mesh.dAdt/mesh.triAreas.^2);
        % dLdt:          [nT        6(v)          3(rgb)    3(phi)    ]
        % mesh.triAreas: [nT                                          ]
        % Lj_1:          [nT                      3(rgb)    3(phi)    ]
        % mesh.dAdt:     [nT        2(xy) 3(vert)                     ]
        % Ki:            [3(phi_col)                        3(phi_row)]
        % dLAdt          [nT        6             3         3         ] 
        % dcdt:          [nT        6             3(rgb)    3(phi)    ]
        dLAdt = (reshape(dLdt,[nT,6,3,3]) - reshape(Lj_1,[nT,1,3,3]).*constantGradients)./mesh.triAreas;
        dcdt = reshape(reshape(dLAdt,[],3)*Ki'/2,nT,6,3,3);
        
        % dcdt:          [nT        6             3(rgb)    3(phi)]
        % Lj_1:          [nT                      3(rgb)    3(phi)]
        % dLdt:          [nT        6             3(rgb)    3(phi)]
        % colors_1:      [nT                      3(rgb)    3(phi)]
        % gradprep:      [nT        6             3(rgb)]
        colors_1 = permute(colors,[1 3 2]);
        % gradPrep_1: [nT 6=(2(xy) x 3(verts)) 3(rgb)]
        gradPrep_1 = sum(dcdt.*reshape(Lj_1,nT,1,3,3) + reshape(colors_1,nT,1,3,3).*dLdt,4);
        % gradPrep: [nT 3(verts) 6=(2(xy) x 3(rgb))] % bundle xy with rgb.
        extra.gradPrep = -permute(reshape(gradPrep_1,nT,2,3,3),[1 3 2 4]);
        
        %% accumulate the per triangle gradient values onto vertices.
        vertGrad = moveTrianglevertValuesToVertexValues(mesh, extra.gradPrep);
        
        % sum over rgb channels
        gradient = sum(vertGrad,3);
        
        % account for mesh boundary
        gradient = slipConditions(mesh, gradient);
    end
end