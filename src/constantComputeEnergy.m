function [energy, colors] = constantComputeEnergy(img, mesh, integral1DNsamples)
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
end