% zero out gradient in direction that violates rectangle boundary.
function grad = polySlipConditions(mesh, grad)
    % slip conds for vertices
    grad.Xgrad(mesh.isXvert,2)=0;
    grad.Xgrad(mesh.isYvert,1)=0;
    
    % slip conds for edges
    grad.edgeGrad(mesh.isXEdge,2,:)=0;
    grad.edgeGrad(mesh.isYEdge,1,:)=0;
    
    % slip conds for faces. noop because no face vertices are on boundary.
end