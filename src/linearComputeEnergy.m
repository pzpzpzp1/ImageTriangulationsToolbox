% colors is (3 vertices) (3 rgb) (nT triangles)
function [energy, colors] = linearComputeEnergy(img, mesh, integral1DNsamples)
    
    X = mesh.X; T = mesh.T; nT = size(T,1);
    % generate sample locations in barycentric coords
    ws = getBarycentricSamplingWeights(integral1DNsamples);
    samplePoints = getSamplePointsFromBarycentricWeights(ws, X, T);
    n = size(ws,1); 
   
    % perform sampling
    sampleVals = double(sampleImage(img, samplePoints)); % n, nT, 3
    
    % compute Lj: nT, 3, 3 
    % triangles, rgb colors, basis funcs
    Lj = squeeze(sum(sampleVals .* reshape(ws,n,1,1,3),1)).*mesh.triAreas/n;
    
    [K, Ki] = loadMassMatrix(1);
    Ki_T = repmat(reshape(Ki,1,3,3),nT,1,1);
    colors = multiprod(permute(Ki_T./(2*mesh.triAreas),[2 3 1]), permute(Lj,[3 2 1]));
    
    % compute actual energy.
    linearColor = permute(reshape(ws * reshape(colors,3,[]), n, 3, nT),[1 3 2]);
    energyPerColorChannel = mesh.triAreas'*squeeze(vecnorm(linearColor - sampleVals,2,1).^2)/n;
    energy = sum(energyPerColorChannel);
end