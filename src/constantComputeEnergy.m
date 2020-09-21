function [extra, energy, colors, gradient] = constantComputeEnergy(img, mesh, integral1DNsamples, salmap)
    extra = {};
    X = mesh.X; T = mesh.T; nT = size(T,1);
    % generate sample locations in barycentric coords
    [ws, interiorInds] = getBarycentricSamplingWeights(integral1DNsamples);
    samplePoints = getSamplePointsFromBarycentricWeights(ws, X, T);
    n = size(ws,1); dA = mesh.triAreas/n;
   
    % perform sampling
    f_triangle = double(sampleImage(img, samplePoints));
    sal_triangle = double(sampleImage(salmap, samplePoints));
    
    % get colors per triangle
    [colors, extra.colorsAlt] = constantGetColorsFromSamples(f_triangle, interiorInds, sal_triangle);
    
    % compute energy integrated over triangles
    extra.perTriangleRGBError = squeeze(sum(sal_triangle.*(f_triangle - reshape(colors,1,nT,3)).^2,1)).*dA;
    energy = sum(extra.perTriangleRGBError,'all');
    
    if nargout >= 4
        %% compute gradient
        % vn is (nT, n, 3, 6). 3 for number of edges. 6 for number of velocity values.
        vndl = sampleVdotN_dl(mesh, integral1DNsamples);

        % generate edge samples of f
        edgeSamplePoints = getTriEdgeSamplePoints(mesh,integral1DNsamples);
        f_tri_edges = double(sampleImage(img, edgeSamplePoints));
        sal_tri_edges = double(sampleImage(salmap, edgeSamplePoints));
        
        % create components needed to compute gradient
        int_svn_dl = squeeze(sum(vndl.*sal_tri_edges, [2, 3]));
        int_s_dA = sum(sal_triangle,1)'.*dA;
        int_fs_dA = colors.*int_s_dA;
        int_fs_dA2 = int_fs_dA.^2;
        int_sfvn_dl = squeeze(sum(vndl.*sal_tri_edges.*reshape(f_tri_edges,nT,integral1DNsamples,3, 1, 3),[2 3]));

        % Build gradient preparation mat. Still vectorized per triangle. Will be re-indexed to lie on vertices.
        % gradPrep: nT(tris per mesh), 3(verts per tri), 2(xy coords), 3(rgb channels)
        A = 2*int_s_dA.*reshape(int_fs_dA,nT,1,3).*int_sfvn_dl;
        B = - reshape(int_fs_dA2,nT,1,3).*int_svn_dl;
        extra.gradPrep = -permute(reshape((A+B)./int_s_dA.^2,nT,2,3,3),[1 3 2 4]);
                
        %% accumulate per triangle gradient values to vertices.
        vertGrad = moveTrianglevertValuesToVertexValues(mesh, extra.gradPrep);
        
        % sum over rgb channels
        gradient = sum(vertGrad,3);

        % account for mesh boundary
        gradient = slipConditions(mesh, gradient);
    end
end