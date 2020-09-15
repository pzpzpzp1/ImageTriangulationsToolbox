% colors is (3 vertices) (3 rgb) (nT triangles)
function [energy, colors, gradient, extra] = linearComputeEnergy(img, mesh, n1D)
    X = mesh.X; T = mesh.T; nT = size(T,1);
    % generate sample locations in barycentric coords
    ws = getBarycentricSamplingWeights(n1D); n = size(ws,1); % number sample points
    samplePoints = getSamplePointsFromBarycentricWeights(ws, X, T);
    n = size(ws,1); 
    
    % perform sampling
    f_triangle = permute(double(sampleImage(img, samplePoints)),[2 1 3]); % nT, n, 3
    
    %{    
        points = reshape(samplePoints,[],2)+[.5 .5];
        cols = reshape(permute(f_triangle,[2 1 3]),[],3);
        figure; hold all; axis equal; 
        image(img);
        scatter(points(:,1),points(:,2),100,cols/255,'filled')
        scatter(points(:,1),points(:,2),'k.')
    %}  
    
    % f_triangle: nT  n     3
    % ws:             n  3
    % compute Lj: nT     3  3 
    % triangles, basis funcs, rgb colors
    Lj = squeeze(sum(reshape(f_triangle,nT,n,1,3) .* reshape(ws,1,n,3,1),2)).*mesh.triAreas/n;
    
    [K, Ki] = loadMassMatrix(1);
    % colors: nT vert rgb
    % Ki     [1   3  3  1 ]
    % Lj_1 = [nT  1  3  3 ]
    % A =    [nT          ]
    extra.colors_fem = squeeze(sum(reshape(Ki,1, 3, 3, 1).*reshape(Lj,nT,1,3,3)./(2*mesh.triAreas), 3));
    
    % use least squares to get color per triangle instead of fem method. fem method doesn't work as well because we use inverse mass matrix which assumes a good numerical integral over the image. turns out uniform samples for integrating on triangles isn't that good.
    A1 = zeros(n,nT,3);    
    A1(:,:,1:2) = samplePoints; A1(:,:,3)=1;
    AtA = squeeze(sum(A1.*reshape(A1,n,nT,1,3),1));
    AtAi = multinv(permute(AtA,[2 3 1]));
    rhsR = squeeze(sum(f_triangle(:,:,1) .* permute(A1,[2 1 3]),2));
    rhsG = squeeze(sum(f_triangle(:,:,2) .* permute(A1,[2 1 3]),2));
    rhsB = squeeze(sum(f_triangle(:,:,3) .* permute(A1,[2 1 3]),2));
    Rcoeffs = squeeze(sum(AtAi .* reshape(rhsR',1,3,nT),2))';
    Gcoeffs = squeeze(sum(AtAi .* reshape(rhsG',1,3,nT),2))';
    Bcoeffs = squeeze(sum(AtAi .* reshape(rhsB',1,3,nT),2))';
    TX = reshape(X(T(:),:),nT,3,2);
    TX1 = TX; TX1(:,:,3)=1;
    rvals = sum(TX1.*reshape(Rcoeffs,nT,1,3),3);
    gvals = sum(TX1.*reshape(Gcoeffs,nT,1,3),3);
    bvals = sum(TX1.*reshape(Bcoeffs,nT,1,3),3);
    colors(:,:,1) = rvals;
    colors(:,:,2) = gvals;
    colors(:,:,3) = bvals;
    
    %{
    figure; hold all; axis equal; 
    image(img);
    renderMeshEdges(mesh,[.5 .5]);
    for i=1:nT
        verts = X(T(i,:),:);
        cent = mean(verts);
        pull = .95;
        for j=1:3
            pt = verts(j,:);
            col = double(squeeze(colors(i,j,:))');
            plotpt = (pull*pt+(1-pull)*cent) + [.5 .5];
            scatter(plotpt(1),plotpt(2),100,col/255,'filled');
            scatter(plotpt(1),plotpt(2),'k.');
        end
    end
    %}
    
    % compute actual energy.
    % ws:         [   n  3  ]
    % colors:     [nT    3 3]
    linearColor = squeeze(sum(reshape(ws,1,n,3,1).*reshape(colors,nT,1,3,3),3));
    energyPerColorChannel = mesh.triAreas'*squeeze(vecnorm(linearColor - f_triangle,2,2).^2)/n;
    energy = sum(energyPerColorChannel);
    
    % todo: optimize. should be slightly redundant with above.
    if nargout >= 3
        %% compute gradient
        
        % constantGradients: [nT 6=(2(xy) 3(verts))]
        constantGradients = mesh.dAdt./mesh.triAreas;
        
        % constantGradients: [nT             6 = (2(xy) 3(verts)) ]
        % ws:                [    n  3(phi)                       ]
        % scratch: dphi_dt(:,:,phind,xyind,vertind) = constantGradients(:,xyind,phind).*ws(:,vertind);
        % dphi_dt: triangles, samples, phi, v
        constantGradients_1 = permute(repmat(reshape(constantGradients,nT,1,2,3),1,n,1,1),[1 2 4 3]); %nT, n, 3 x 2.
        ws_1 = zeros(nT,n,2,2,3);
        ws_1(:,:,1,1,:) = repmat(reshape(ws,1,n,3),nT,1,1);
        ws_1(:,:,2,2,:) = repmat(reshape(ws,1,n,3),nT,1,1);
        ws_2 = reshape(ws_1,nT,n,2,6); % nT n 2 6
        dphi_dt = squeeze(sum(constantGradients_1.*reshape(ws_2,nT,n,1,2,6),4));

        % vn is (nT, n, 3, 6). 3 for number of edges. 6 for number of velocity values.
        vndl = sampleVdotN_dl(mesh, n1D);

        % generate edge samples of f.
        % edgeSamplePoints: [nT, n, 3(edges), 2(xy)]
        edgeSamplePoints = getEdgeSamplePoints(mesh,n1D);
        % f_tri_edges:      [nT, n, 3(edges), 3(rgb)]
        f_tri_edges = double(sampleImage(img, edgeSamplePoints));

%         figure; hold all; axis equal; image(img);
%         pts = reshape(edgeSamplePoints,[],2);
%         cols = double(reshape(f_tri_edges,[],3))/255;
%         scatter(pts(:,1),pts(:,2),200,cols,'filled')
%         scatter(pts(:,1),pts(:,2),'k.')
        
        % generate edge samples of phi
        % phi_edges_oneTri: (samples per edge) (3 edges) (3 phis)
        edgeWs = linspace(0,1,n1D);
        phi_edges_oneTri = zeros(n1D,3,3);
        phi_edges_oneTri(:,1,1) = fliplr(edgeWs);
        phi_edges_oneTri(:,2,2) = fliplr(edgeWs);
        phi_edges_oneTri(:,3,3) = fliplr(edgeWs);
        phi_edges_oneTri(:,3,1) = edgeWs;
        phi_edges_oneTri(:,1,2) = edgeWs;
        phi_edges_oneTri(:,2,3) = edgeWs;
        
        % phi_edges_oneTri: [      n1D    3(edges)                 3(phis)]
        % vndl:             [nT    n1D    3(edges)   6(v)                 ]
        % f_tri_edges:     [nT    n1D    3(edges)          3(rgb)        ]
        % int_fvn_dl: (nT) (v motions) (rgb) (phi)
        int_fvnphi_dl = squeeze(sum(reshape(phi_edges_oneTri,1,n1D,3,1,1,3).*vndl.*reshape(f_tri_edges,nT,n1D,3,1,3,1),[2 3]));
        
        % build: dLdt
        % f_triangle_1:  [nT n   3(rgb)         ]
        % dphi_dt_1:     [nT n 6         3(phi) ]
        % int_fvnphi_dl: [nT   6 3(rgb)  3(phi) ]
        % dLdt:          [nT   6 3(rgb)  3(phi) ]
        dphi_dt_1 = reshape(permute(dphi_dt,[1 2 4 3]),nT,n,6,1,3);
        f_triangle_1 = reshape(f_triangle,nT,n,1,3,1);
        dLdt = squeeze(sum(dphi_dt_1.*f_triangle_1.*(mesh.triAreas/n),2)) + int_fvnphi_dl;
        
        % build: dcdt
        % dLdt:          [nT        6(v)          3(rgb)    3(phi)    ]
        % mesh.triAreas: [nT]
        % Lj_1:          [nT                      3(rgb)    3(phi)    ]
        % mesh.dAdt:     [nT        2(xy) 3(vert)                     ]
        % Ki:            [3(phi_col)                        3(phi_row)]
        % dLAdt          [nT        6             3         3         ]
        % dcdt pseudocode = Ki/2 * (dLdt/mesh.triAreas - Lj.*mesh.dAdt/mesh.triAreas.^2);
        % dcdt:          [nT        6             3(rgb)    3(phi)    ]
        Lj_1 = permute(Lj,[1 3 2]);
        dLAdt = reshape(dLdt,[nT,6,3,3])./mesh.triAreas - reshape(Lj_1,[nT,1,3,3]).*mesh.dAdt./mesh.triAreas.^2;
        dcdt = reshape(reshape(dLAdt,[],3)*Ki'/2,nT,6,3,3);
        
        %% grad = d/dt int_T f^2 when f is a constant gradient. verified!
        int_f2vn_dl = squeeze(sum(vndl.*reshape(f_tri_edges.^2,nT,n1D,3,1,3,1),[2 3]));
        %% dgdt = 0 when f is constant gradient
        % colors:      [nT             3(verts)    3(rgb)
        % dphi_dt_1:   [nT   n    (2 x 3(verts))            3(phi)
        % dcdt:        [nT        (2 x 3(verts))   3(rgb)   3(phi)    
        % phi:         [     n                              3(phi)    
        dphi_dt_1 = permute(dphi_dt,[1 2 4 3]);
        dgdt = sum(reshape(colors,nT,1,1,3,3).*reshape(dphi_dt,nT,n,2,3,1,3) + reshape(dcdt,nT,1,2,3,3,3).*reshape(ws,1,n,1,1,1,3),6);
        dgdt = reshape(dgdt,nT,n,6,3);
        
        % dcdt:          [nT        6             3(rgb)    3(phi)]
        % Lj:            [nT                      3(rgb)    3(phi)]
        % dLdt:          [nT        6             3(rgb)    3(phi)]
        % colors_1:      [nT                      3(rgb)    3(phi)]
        % gradprep:      [nT        6             3(rgb)]
        colors_1 = permute(colors,[1 3 2]);
        % gradPrep_1: [nT 6=(2(xy) x 3(verts)) 3(rgb)]
        gradPrep_1 = sum(dcdt.*reshape(permute(Lj,[1,3,2]),nT,1,3,3) + reshape(colors_1,nT,1,3,3).*dLdt,4);
        gradPrep_1 = int_f2vn_dl;
        % gradPrep: [nT 3(verts) 6=(2(xy) x 3(rgb))]
        gradPrep = -reshape(permute(reshape(gradPrep_1,nT,2,3,3),[1 3 2 4]),nT,3,6); % bundle xy with rgb.
        
        
        % accumulate per triangle gradient values to vertices.
        nX = size(X,1); 
        vertGrad = zeros(nX,6);
        for i=1:6
            vertGrad(:,i) = accumarray(mesh.T(:), reshape(gradPrep(:,:,i),[],1));
        end
        vertGrad = reshape(vertGrad,nX,2,3);

        % sum over rgb channels
        gradient = sum(vertGrad,3);
        assert(norm(reshape(vertGrad(:,:,1)-vertGrad(:,:,2),[],2))==0);
        
        % account for mesh boundary
        gradient = slipConditions(mesh, gradient);
    end
end