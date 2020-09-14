function [energy, colors, gradient, extra] = constantComputeEnergy(img, mesh, integral1DNsamples)
    extra = {};
    X = mesh.X; T = mesh.T; nT = size(T,1);
    % generate sample locations in barycentric coords
    ws = getBarycentricSamplingWeights(integral1DNsamples);
    samplePoints = getSamplePointsFromBarycentricWeights(ws, X, T);
    n = size(ws,1); 
   
    % perform sampling
    sampleVals = sampleImage(img, samplePoints);
    
    % compute mean per triangle
    colorsDouble = (squeeze(mean(sampleVals,1)));
    colors = uint8(ceil(squeeze(mean(sampleVals,1))));
    
    % compute energy integrated over triangles
    perTriangleRGBError = squeeze(sum((double(sampleVals) - reshape(colorsDouble,1,nT,3)).^2,1));
    energy = sum(mesh.triAreas' * perTriangleRGBError/n);
    
    % Todo: optimize gradient computation. Should be slightly redundant with energy part.
    if nargout >= 3
        %% compute gradient
        X = mesh.X; T = mesh.T; edges = mesh.edges;
        nX = size(X,1); nT = size(T,1); nE = size(edges,1);

        % generate area samples of f
        ws = getBarycentricSamplingWeights(integral1DNsamples); n = size(ws,1); % number sample points
        samplePoints = getSamplePointsFromBarycentricWeights(ws, X, T); % n x nT x 2
        f_triangle = sampleImage(img, samplePoints); % n x nT x 3

        % vn is (nT, n, 3, 6). 3 for number of edges. 6 for number of velocity values.
        vndl = sampleVdotN_dl(mesh, integral1DNsamples);

        % generate edge samples of f
        edgeSamplePoints = getEdgeSamplePoints(mesh,integral1DNsamples);
        f_tri_edges = sampleImage(img, edgeSamplePoints);

        % create components needed to compute gradient
        int_f_dA = squeeze(sum(f_triangle,1)).*mesh.triAreas/n;
        int_f_dA2 = int_f_dA.^2;
        % int_vn_dl = squeeze(sum(vndl, [1, 3]));
        int_vn_dl = mesh.dAdt;
        int_fvn_dl = squeeze(sum(vndl.*reshape(double(f_tri_edges),nT,integral1DNsamples,3, 1, 3),[2 3]));

        % Build gradient preparation mat. Still vectorized per triangle. Will be re-indexed to lie on vertices.
        % gradPrep: nT(tris per mesh), 3(verts per tri), 2(xy coords), 3(rgb channels)
        A = 2*mesh.triAreas.*repmat(reshape(int_f_dA,nT,1,3),1,6).*int_fvn_dl;
        B = - repmat(reshape(int_f_dA2,nT,1,3),1,6).*repmat(int_vn_dl,1,1,3);
        gradPrep = -reshape(permute(reshape((A+B)./mesh.triAreas.^2,nT,2,3,3),[1 3 2 4]),nT,3,6); 

        % accumulate per triangle gradient values to vertices.
        vertGrad = zeros(nX,6);
        for i=1:6
            vertGrad(:,i) = accumarray(mesh.T(:), reshape(gradPrep(:,:,i),[],1));
        end
        vertGrad = reshape(vertGrad,nX,2,3);

        % sum over rgb channels
        gradient = sum(vertGrad,3);

        % account for mesh boundary
        gradient = slipConditions(mesh, gradient);
    end
end