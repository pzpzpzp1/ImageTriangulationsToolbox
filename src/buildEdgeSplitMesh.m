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
    ip1 = mod(JJ(perm)+1,3); ip1(ip1==0)=3;
    ip2 = mod(JJ(perm)+2,3); ip2(ip2==0)=3;
    v2 = v123(sub2ind(size(v123),[1:nE]',ip1));
    v3 = v123(sub2ind(size(v123),[1:nE]',ip2));
    [II,JJ] = find(~squeeze(sum(reshape(v234,nE,1,3) == v23,2)));
    [~,perm] = sort(II);
    v4 = v234(sub2ind(size(v234),[1:nE]',JJ(perm)));
    
    % create new vertices at the center of each edge;
    v5 = [(nX+1):(nX+nE)]';
    Xnew = [X; (X(mesh.edges(:,1),:)+X(mesh.edges(:,2),:))/2];
    % build triangulation
    Tnew1 = [v1 v2 v5; v1 v5 v3; v5 v2 v4; v3 v5 v4];
    Tnew = reshape(permute(reshape(Tnew1,nE,4,3),[2 1 3]),[],3);
    
    % reorient inverted triangles
    flippedTris = getTriangleAreas(Xnew, Tnew) < 0;
    Tnew(flippedTris,:) = Tnew(flippedTris,[1 3 2]);
    
    % triangle soup
%     figure; hold all; axis equal;
%     patch('vertices',Xnew,'faces',Tnew(1:4:end,:),'facecolor','green','facealpha',.1);
%     patch('vertices',Xnew,'faces',Tnew(2:4:end,:),'facecolor','red','facealpha',.1);
%     patch('vertices',Xnew,'faces',Tnew(3:4:end,:),'facecolor','blue','facealpha',.1);
%     patch('vertices',Xnew,'faces',Tnew(4:4:end,:),'facecolor','yellow','facealpha',.1);
    
    esMesh = MeshFromXT(Xnew, Tnew, 1);
end




