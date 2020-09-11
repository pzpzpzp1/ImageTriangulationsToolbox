% zero out gradient in direction that violates rectangle boundary.
function grad = slipConditions(mesh, grad)
    grad(mesh.isXvert,2)=0;
    grad(mesh.isYvert,1)=0;
end