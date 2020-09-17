function esMesh = buildEdgeSplitMesh(mesh)
    nE = mesh.nE; nT = mesh.nT; nX = mesh.nX;
    T = mesh.T; X = mesh.X;
    
    % each edge is made of verts v23. opposing vertices are v1 and v3
    v123 = T(mesh.edges2triangles(:,1),:);
    v234 = T(mesh.edges2triangles(:,2),:);
    v23 = mesh.edges;
    [II,JJ] = find(~squeeze(sum(reshape(v123,nE,1,3) == v23,2)));
    [~,perm] = sort(II);
    v1 = v123(sub2ind(size(v123),[1:nE]',JJ(perm)));
    [II,JJ] = find(~squeeze(sum(reshape(v234,nE,1,3) == v23,2)));
    [~,perm] = sort(II);
    v4 = v234(sub2ind(size(v234),[1:nE]',JJ(perm)));
    v2 = v23(:,1);
    v3 = v23(:,2);
    
    % create new vertices at the center of each edge;
    v5 = [(nX+1):(nX+nE)]';
    Xnew = [X; (X(mesh.edges(:,1),:)+X(mesh.edges(:,2),:))/2];
    % build triangulation
    Tnew1 = [v1 v2 v5; v5 v2 v4; v3 v5 v4; v1 v5 v3];
    Tnew = reshape(permute(reshape(Tnew1,nE,4,3),[2 1 3]),[],3);
    
    % reorient inverted triangles
    v1 = Xnew(Tnew(:,1),:);  v2 = Xnew(Tnew(:,2),:);  v3 = Xnew(Tnew(:,3),:);
    e12 = [v1-v2]; e12(1,3) = 0; e23 = [v2-v3]; e23(1,3) = 0;
    flippedTris = (cross(e12,e23)/2)*[0 0 1]' < 0;
    Tnew(flippedTris,:) = Tnew(flippedTris,[1 3 2]);
    
    % triangle soup
%     figure; hold all; axis equal;
%     patch('vertices',Xnew,'faces',Tnew(1:4:end,:),'facecolor','green','facealpha',.1);
%     patch('vertices',Xnew,'faces',Tnew(2:4:end,:),'facecolor','red','facealpha',.1);
%     patch('vertices',Xnew,'faces',Tnew(3:4:end,:),'facecolor','blue','facealpha',.1);
%     patch('vertices',Xnew,'faces',Tnew(4:4:end,:),'facecolor','yellow','facealpha',.1);
    
    esMesh = MeshFromXT(Xnew, Tnew, 1);
end




