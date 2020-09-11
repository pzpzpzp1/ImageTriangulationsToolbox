function mesh = MeshFromXT(X,T)
    mesh.X = X;
    mesh.T = T;
    mesh.nX = size(X,1);
    mesh.nT = size(T,1);
    
    v1 = X(T(:,1),:);  v2 = X(T(:,2),:);  v3 = X(T(:,3),:);
    e12 = [v1-v2]; e12(1,3) = 0; e23 = [v2-v3]; e23(1,3) = 0;
    mesh.triAreas = (cross(e12,e23)/2)*[0 0 1]';
    
    allEdges = [T(:,[1 2]); T(:,[2 3]); T(:,[3 1])];
    [mesh.edges, ia, ic] = unique(sort(allEdges,2),'rows');
    
    mesh.triangles2edges = reshape(ic,size(T,1),3);
    
    e1 = mesh.edges(mesh.triangles2edges(:,1),:);
    e2 = mesh.edges(mesh.triangles2edges(:,2),:);
    e3 = mesh.edges(mesh.triangles2edges(:,3),:);
    
    e123 = reshape(mesh.edges(mesh.triangles2edges(:),:), mesh.nT, 3, 2);
    tang123 = reshape(X(e123(:,:,2),:)-X(e123(:,:,1),:),mesh.nT,3,2);
    mesh.triangleEdgeNormals = zeros(mesh.nT,3,2);
    mesh.triangleEdgeNormals(:,:,1) = - tang123(:,:,2);
    mesh.triangleEdgeNormals(:,:,2) = tang123(:,:,1);
    mesh.triangleEdgeLengths = (vecnorm(tang123,2,3));
    
    
    
end