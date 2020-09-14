function mesh = MeshFromXT(X,T)
    mesh.X = X;
    mesh.T = T;
    mesh.nX = size(X,1);
    mesh.nT = size(T,1);
    
    v1 = X(T(:,1),:);  v2 = X(T(:,2),:);  v3 = X(T(:,3),:);
    e12 = [v1-v2]; e12(1,3) = 0; e23 = [v2-v3]; e23(1,3) = 0;
    mesh.triAreas = (cross(e12,e23)/2)*[0 0 1]';
    assert(all(mesh.triAreas>0));
    
    allEdges = [T(:,[1 2]); T(:,[2 3]); T(:,[3 1])];
    [mesh.edges, ia, ic] = unique(sort(allEdges,2),'rows');
    
    %[~, inds] = ismember(reshape(sort(permute(reshape(mesh.T(:,[1 2 2 3 3 1]),[],2,3),[1 3 2]),3),mesh.nT*3,2),  sort(mesh.edges,2), 'rows');
    % mesh.triangles2edges = reshape(inds,mesh.nT,3);
    mesh.triangles2edges = reshape(ic,size(T,1),3);
    
    e1 = v2-v1;
    e2 = v3-v2;
    e3 = v1-v3;
    mesh.triangleEdgeNormals = zeros(mesh.nT,3,2);
    mesh.triangleEdgeNormals(:,1,1) = - e1(:,2);
    mesh.triangleEdgeNormals(:,1,2) = e1(:,1);
    mesh.triangleEdgeNormals(:,2,1) = - e2(:,2);
    mesh.triangleEdgeNormals(:,2,2) = e2(:,1);
    mesh.triangleEdgeNormals(:,3,1) = - e3(:,2);
    mesh.triangleEdgeNormals(:,3,2) = e3(:,1);
    mesh.triangleEdgeNormals = -mesh.triangleEdgeNormals ./ vecnorm(mesh.triangleEdgeNormals,2,3);
    mesh.triangleEdgeLengths = [vecnorm(e1,2,2) vecnorm(e2,2,2) vecnorm(e3,2,2)];
    
%     figure; hold all;axis equal;
%     for i=1:mesh.nT
%         triv = mesh.X(mesh.T(i,:),:);
%         triv = triv*.9 + .1*mean(triv);
%         ecs = (triv + circshift(triv,-1)) /2;
%         patch('vertices',triv,'faces',[1 2 3],'facecolor','green','facealpha',.2)
%         en = squeeze(mesh.triangleEdgeNormals(i,:,:));
%         quiver(ecs (:,1),ecs (:,2),en(:,1),en(:,2),'r');
%         pause
%     end
    
    %% boundary handling. Specific to rectangular boundary mesh!
    isBoundaryEdge = accumarray(mesh.triangles2edges(:),ones(3*size(T,1),1))==1;
    edgetangs = X(mesh.edges(:,1),:)-X(mesh.edges(:,2),:);
    edgetangs = edgetangs./vecnorm(edgetangs,2,2);
    isXEdge = abs(edgetangs*[0 1]')<.000001;
    isYEdge = abs(edgetangs*[1 0]')<.000001;
    XvertInds = unique(mesh.edges(isBoundaryEdge & isXEdge,:));
    YvertInds = unique(mesh.edges(isBoundaryEdge & isYEdge,:));
    mesh.isXvert = false(size(X,1),1); mesh.isXvert(XvertInds) = true;
    mesh.isYvert = false(size(X,1),1); mesh.isYvert(YvertInds) = true;
    mesh.isInterior = ~mesh.isXvert & ~mesh.isYvert;
    
    %% compute dA_dt
    TX = permute(reshape(X(T',:),3,[],2),[2 3 1]);
    e123 = TX - circshift(TX,-1,3);
    e312 = circshift(e123,-1,3); % get part of e123 orthogonal to e312
    altitudes = e123 - sum(e123.*e312,2)./vecnorm(e312,2,2).^2.*e312;
    altitudes = altitudes./vecnorm(altitudes,2,2);
    mesh.dAdt = reshape(vecnorm(e312,2,2).*altitudes/2,[],6); %nT, 6=(2(xy) x 3(verts)).
    
    %{
    dirs = squeeze(altitudes(end,:,:))';
    XT = X(T(end,:),:);
    figure; hold all; axis equal;
    plot(XT(:,1), XT(:,2));
    quiver(XT(:,1), XT(:,2), dirs(:,1), dirs(:,2));
    %}
    
    
    
end