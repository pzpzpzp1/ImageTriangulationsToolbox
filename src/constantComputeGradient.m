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
    f_edges = zeros(nT,integral1DNsamples,3,2);
    edgeSamplePoints = getEdgeSamplePoints(mesh,integral1DNsamples);
    f_tri_edges = sampleImage(img, edgeSamplePoints);
    
    
    
    
    int_f_dA = squeeze(sum(f_triangle,1)).*mesh.triAreas/n;
    int_vn_dl = squeeze(sum(vndl, [1, 3]));
    
    
end