function [energyPerTriangle] = polyGetAreaEnergyCore(mesh, integral1DNsamples, factormap)
    
    ws = getBarycentricSamplingWeights(integral1DNsamples);
    [samplePoints, ~, geodets] = mesh.triMaps(ws); % world space points, and jacobian det
    factorSamples = double(sampleImage(factormap, samplePoints)); % factor map values
    energyPerTriangle = -dot(log(geodets), factorSamples, 1);
    
end