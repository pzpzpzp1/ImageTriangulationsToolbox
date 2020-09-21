function [energy, gradient] = getAreaEnergy(mesh)
    nT = size(mesh.T,1);
    energy = sum(-log(mesh.triAreas));
    gradprep = permute(-reshape(mesh.dAdt,nT,2,3)./mesh.triAreas, [1 3 2]);
    vertgrad = moveTrianglevertValuesToVertexValues(mesh, gradprep);
    gradient = slipConditions(mesh, vertgrad);
end