function gradient = constantComputeGradient(img, mesh, integral1DNsamples)
    X = mesh.X; T = mesh.T; edges = mesh.edges;
    nX = size(X,1); nT = size(T,1); nE = size(edges,1);
    
    % generate area samples of f
    ws = getBarycentricSamplingWeights(integral1DNsamples); n = size(ws,1); % number sample points
    samplePoints = getSamplePointsFromBarycentricWeights(ws, X, T); % n x nT x 2
    f_triangle = sampleImage(img, samplePoints); % n x nT x 3
    
    % vn is (n, nT, 3, 6). 3 for number of edges. 6 for number of velocity values.
    vndl = sampleVdotN_dl(mesh, integral1DNsamples);
    
    % generate edge samples of f
    edgeSamplePoints = getEdgeSamplePoints(mesh,integral1DNsamples);
    f_tri_edges = sampleImage(img, edgeSamplePoints);
    
    % create components needed to compute gradient
    int_f_dA = squeeze(sum(f_triangle,1)).*mesh.triAreas/n;
    int_f_dA2 = int_f_dA.^2;
    % int_vn_dl = squeeze(sum(vndl, [1, 3]));
    int_vn_dl = reshape(mesh.dAdt,[],6);
    int_fvn_dl = squeeze(sum(permute(vndl,[2 1 3 4]).*reshape(double(f_tri_edges),nT,integral1DNsamples,3, 1, 3),[2 3]));
    
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