function esMesh = buildTriSplitMesh(mesh)
    if nargin == 0
        [X,T] = initialGridMesh(10, 10, 5, 1);
        mesh = MeshFromXT(X,T);
    end
    X = mesh.X; T = mesh.T;
    
    newX = (X(mesh.edges(:,1),:)+X(mesh.edges(:,2),:))/2;
    Xnew = [X; newX];
    newvinds = (mesh.nX+1):size(Xnew,1);
    
    v1 = mesh.T(:,1);     v2 = mesh.T(:,2);    v3 = mesh.T(:,3);
    v4 = newvinds(mesh.triangles2edges(:,1))';
    v5 = newvinds(mesh.triangles2edges(:,2))';
    v6 = newvinds(mesh.triangles2edges(:,3))';
    
    Tnew = [v1 v4 v6; v4 v2 v5; v6 v5 v3; v6 v4 v5];
    Tnew = reshape(permute(reshape(Tnew,mesh.nT,4,3),[2 1 3]),[],3);
    
    esMesh = MeshFromXT(Xnew,Tnew);
end