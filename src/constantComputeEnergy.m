function [extra, energy, colors, gradient] = constantComputeEnergy(img, mesh, integral1DNsamples)
    extra = {};
    X = mesh.X; T = mesh.T; nT = size(T,1);
    % generate sample locations in barycentric coords
    [ws, interiorInds] = getBarycentricSamplingWeights(integral1DNsamples);
    samplePoints = getSamplePointsFromBarycentricWeights(ws, X, T);
    n = size(ws,1); dA = mesh.triAreas/n;
   
    % perform sampling
    f_triangle = double(sampleImage(img, samplePoints));
    
    % get colors per triangle
    [colors, extra.colorsAlt] = constantGetColorsFromSamples(f_triangle, interiorInds);
    
    % compute energy integrated over triangles
    extra.perTriangleRGBError = squeeze(sum((f_triangle - reshape(colors,1,nT,3)).^2,1)).*dA;
    energy = sum(extra.perTriangleRGBError,'all');
    
    if nargout >= 4
        %% compute gradient
        % vn is (nT, n, 3, 6). 3 for number of edges. 6 for number of velocity values.
        vndl = sampleVdotN_dl(mesh, integral1DNsamples);

        % generate edge samples of f
        edgeSamplePoints = getEdgeSamplePoints(mesh,integral1DNsamples);
        f_tri_edges = double(sampleImage(img, edgeSamplePoints));

        % create components needed to compute gradient
        int_f_dA = colors.*mesh.triAreas;
        int_f_dA2 = int_f_dA.^2;
        int_fvn_dl = squeeze(sum(vndl.*reshape(f_tri_edges,nT,integral1DNsamples,3, 1, 3),[2 3]));

        % Build gradient preparation mat. Still vectorized per triangle. Will be re-indexed to lie on vertices.
        % gradPrep: nT(tris per mesh), 3(verts per tri), 2(xy coords), 3(rgb channels)
        A = 2*mesh.triAreas.*reshape(int_f_dA,nT,1,3).*int_fvn_dl;
        B = - reshape(int_f_dA2,nT,1,3).*mesh.dAdt;
        extra.gradPrep = -permute(reshape((A+B)./mesh.triAreas.^2,nT,2,3,3),[1 3 2 4]);
        
        %% accumulate per triangle gradient values to vertices.
        vertGrad = moveTrianglevertValuesToVertexValues(mesh, extra.gradPrep);
        
        % sum over rgb channels
        gradient = sum(vertGrad,3);

        % account for mesh boundary
        gradient = slipConditions(mesh, gradient);
    end
end