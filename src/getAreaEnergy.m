function [energy, gradient] = getAreaEnergy(mesh, salmap)
    nT = size(mesh.T,1);
    energy = sum(-log(mesh.triAreas));
    gradprep = permute(-reshape(mesh.dAdt,nT,2,3)./mesh.triAreas, [1 3 2]);
    
    if numel(salmap)~=0
        % weight triangle graident by saliency map
        [ws, interiorInds] = getBarycentricSamplingWeights(5);
        samplePoints = getSamplePointsFromBarycentricWeights(ws, mesh.X, mesh.T);
        sal_triangle = double(sampleImage(salmap, samplePoints));
        salgradprep = gradprep ./ mean(sal_triangle)';
    else
        salgradprep = gradprep;
    end
    
    vertgrad = moveTrianglevertValuesToVertexValues(mesh, salgradprep);
    gradient = slipConditions(mesh, vertgrad);
end