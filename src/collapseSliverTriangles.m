% given list of presumably sliver triangle indices, collapse them. 
% find longest edge and project opposite vertex onto it. The edge that's projected on gets split and its opposing triangle is also split into two.
% boundary case and interior cases handled separately both are vectorized... although it may be unnecessary.
function [esmesh, reducedInds] =  collapseSliverTriangles(mesh, badTriInds)
    if nargin==0
        load temp.mat;
        T = mesh.T; X = mesh.X;
        mesh = MeshFromXT(X,T);
        % save('temp.mat','mesh','badTriInds');
    end
    T = mesh.T; X = mesh.X;
    Xnew = X;
    
    TEL = mesh.triangleEdgeLengths(badTriInds,:)';
    TE = mesh.triangles2edges(badTriInds,:)';
    [~,maxEdgeIndLinear] = max(TEL,[],'linear');
    longEdgeInd = TE(maxEdgeIndLinear)';
    
    %% get subset that don't conflict
    isCoveredVertices = false(mesh.nX,1);
    keepEdgeInd = false(numel(longEdgeInd),1);
    for i=1:numel(longEdgeInd)
        eind = longEdgeInd(i);
        ts = mesh.T(mesh.edges2triangles(eind,:),:);
        if ~any(isCoveredVertices(ts))
            isCoveredVertices(ts) = true;
            keepEdgeInd(i) = true;
        end
    end
    badTriInds = badTriInds(keepEdgeInd);
    longEdgeInd = longEdgeInd(keepEdgeInd);

    %% boundary edge case
    inds = find(mesh.isBoundaryEdge(longEdgeInd));
    boundaryLongEdge = longEdgeInd(inds);
    boundaryBadTri = badTriInds(inds);
    ne = numel(boundaryLongEdge);
    v123 = T(mesh.edges2triangles(boundaryLongEdge,1),:);
    v23 = mesh.edges(boundaryLongEdge,:);
    [II,JJ] = find(~permute(sum(reshape(v123, numel(boundaryLongEdge) ,1,3) == v23,2),[1 3 2]));
    [~,perm] = sort(II);
    v1 = v123(sub2ind(size(v123),[1:ne]',JJ(perm)));
    triangles2remove = boundaryBadTri(:);
    v2 = v23(:,1);
    Xnew(v1(mesh.isXEdge(boundaryLongEdge)),2)=Xnew(v2(mesh.isXEdge(boundaryLongEdge)),2);
    Xnew(v1(mesh.isYEdge(boundaryLongEdge)),1)=Xnew(v2(mesh.isYEdge(boundaryLongEdge)),1);
    
    %% interior edge case. also, remove boundary triangles whose largest edge isn't boundary.
    inds = find(~mesh.isBoundaryEdge(longEdgeInd) & ~mesh.isBoundaryTriangle(badTriInds));
    boundaryLongEdge = longEdgeInd(inds);
    boundaryBadTri = badTriInds(inds);
    ne = numel(boundaryLongEdge);
    v23 = mesh.edges(boundaryLongEdge,:);
    v123 = T(mesh.edges2triangles(boundaryLongEdge,1),:);
    v234 = T(mesh.edges2triangles(boundaryLongEdge,2),:);
    [II,JJ] = find(~permute(sum(reshape(v123, numel(boundaryLongEdge) ,1,3) == v23,2),[1 3 2]));
    [~,perm] = sort(II);
    v1 = v123(sub2ind(size(v123),[1:ne]',JJ(perm)));
    ip1 = mod(JJ(perm)+1,3); ip1(ip1==0)=3;
    ip2 = mod(JJ(perm)+2,3); ip2(ip2==0)=3;
    v2 = v123(sub2ind(size(v123),[1:ne ]',ip1));
    v3 = v123(sub2ind(size(v123),[1:ne ]',ip2));
    [II,JJ] = find(~permute(sum(reshape(v234, numel(boundaryLongEdge) ,1,3) == v23,2),[1 3 2]));
    [~,perm] = sort(II);
    v4 = v234(sub2ind(size(v234),[1:ne]',JJ(perm)));
    triangles2remove = unique([triangles2remove; reshape(unique(mesh.edges2triangles(boundaryLongEdge,:)),[],1)]);
    Tnew = [T; v1 v4 v3; v1 v2 v4];
    v1inds = any((T(boundaryBadTri,:) - v1)==0,2);
    v4inds = any((T(boundaryBadTri,:) - v4)==0,2);
    Xnew(v4(v4inds),:) = projectPointOntoLine(Xnew(v4(v4inds),:), Xnew(v2(v4inds),:), Xnew(v3(v4inds),:));
    Xnew(v1(v1inds),:) = projectPointOntoLine(Xnew(v1(v1inds),:), Xnew(v2(v1inds),:), Xnew(v3(v1inds),:));
    
    %% build final mesh
    Tnew(triangles2remove,:)=[];
    
    %% check for inversions from collapse
    % the epsilon area buffer is needed becasue in trisplit step, areas are cut even further. sometimes that fails at floating precision and turns non inverted triangles into inverted ones.
    flipped = find(getTriangleAreas(Xnew,Tnew)<=1e-9); 
    reducedInds = [];
    if numel(flipped)==0;
        esmesh = MeshFromXT(Xnew, Tnew);
    else
        esmesh = mesh;
    
        reducedInds = badTriInds(~any(ismember(T(badTriInds,:), Tnew(flipped,:)),2));
    end
    
    %{
    % before
    figure; image(img)
    ax1=subplot(1,2,1);
    axis equal; hold all;
    patch('faces',T,'vertices',X,'facecolor','green','facealpha',.1)
    patch('faces',T(badTriInds,:),'vertices',X,'facecolor','none','facealpha',1,'edgecolor','red','linewidth',2)
    
    % after
    ax2=subplot(1,2,2);
    axis equal; hold all;
    patch('faces',Tnew,'vertices',Xnew,'facecolor','green','facealpha',.1)
    patch('faces',Tnew(flipped,:),'vertices',Xnew,'facecolor','cyan','edgecolor','cyan')
    
    linkprop([ax1 ax2],{'xlim','ylim','zlim'});
    %}
    
end