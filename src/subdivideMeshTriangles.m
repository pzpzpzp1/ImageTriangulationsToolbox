function [Xnew, Tnew] = subdivideMeshTriangles(mesh, triInds)
    if nargin == 0
        [X,T] = initialGridMesh(10, 10, 7, 1);
        mesh = MeshFromXT(X,T);
        triInds = randsample(1:mesh.nT, 10)';
    end
    X = mesh.X; T = mesh.T; nT = size(T,1); nX = size(X,1);
    
    alsosplit = true;
    while any(alsosplit)
        splitEdgeInds = unique(mesh.triangles2edges(triInds,:));
        splitedgecount = sum(ismember(mesh.triangles2edges, splitEdgeInds),2);
        alsosplit = splitedgecount == 2;
        triInds = unique([triInds; find(alsosplit)]);
    end
    
    noSplitTris = splitedgecount==0;
    oneSplitTris = splitedgecount==1;
    fullSplitTris = splitedgecount==3;
    
    e2e = zeros(mesh.nE,1);
    edges2split = unique(mesh.triangles2edges(triInds,:));
    e2e(edges2split) = 1:numel(edges2split);
    
    newX = (X(mesh.edges(edges2split,1),:)+X(mesh.edges(edges2split,2),:))/2;
    Xnew = [X; newX];
    newvinds = (mesh.nX+1):size(Xnew,1);
    
    T0new = [T(noSplitTris,:);];
    
    v1 = mesh.T(fullSplitTris,1);     v2 = mesh.T(fullSplitTris,2);    v3 = mesh.T(fullSplitTris,3);
    v4 = newvinds(e2e(mesh.triangles2edges(fullSplitTris,1)))';
    v5 = newvinds(e2e(mesh.triangles2edges(fullSplitTris,2)))';
    v6 = newvinds(e2e(mesh.triangles2edges(fullSplitTris,3)))';
    T3new = [v1 v4 v6; v4 v2 v5; v6 v5 v3; v6 v4 v5];
    
    e123 = mesh.triangles2edges(oneSplitTris,:);
    [~,d] = max(ismember(e123, edges2split)');
    oneSplitTriInds = find(oneSplitTris);
    e1inds = (oneSplitTriInds(d==1));
    e2inds = (oneSplitTriInds(d==2));
    e3inds = (oneSplitTriInds(d==3));
    
    ve1 = newvinds(e2e(mesh.triangles2edges(e1inds,1)))';
    ve2 = newvinds(e2e(mesh.triangles2edges(e2inds,2)))';
    ve3 = newvinds(e2e(mesh.triangles2edges(e3inds,3)))';
    T1new1 = [mesh.T(e1inds,1) ve1 mesh.T(e1inds,3); mesh.T(e1inds,3) ve1 mesh.T(e1inds,2);];
    T1new2 = [mesh.T(e2inds,1) mesh.T(e2inds,2) ve2; mesh.T(e2inds,1) ve2 mesh.T(e2inds,3);];
    T1new3 = [mesh.T(e3inds,1) mesh.T(e3inds,2) ve3; ve3 mesh.T(e3inds,2) mesh.T(e3inds,3);];
    
    Tnew = [T0new; T3new; T1new1; T1new2; T1new3;];
    
    %{
    figure; hold all; axis equal; 
    patch('faces',Tnew,'vertices',Xnew,'facecolor','blue','facealpha',.5);    
    patch('faces',T0new,'vertices',Xnew,'facecolor','green','facealpha',.5);
    patch('faces',T3new,'vertices',Xnew,'facecolor','red','facealpha',.5);
    patch('faces',T1new1,'vertices',Xnew,'facecolor','blue','facealpha',.5);
    patch('faces',T1new2,'vertices',Xnew,'facecolor','blue','facealpha',.5);
    patch('faces',T1new3,'vertices',Xnew,'facecolor','blue','facealpha',.5);
    assert(all(getTriangleAreas(Xnew,Tnew) > 0))
    %}
end