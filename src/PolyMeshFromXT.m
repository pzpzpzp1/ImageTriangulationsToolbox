function mesh = PolyMeshFromXT(X,T,edges,edgeX,faceX,isTriangleSoup,gdeg)
    if nargin==0
        [X,T] = initialHexLatticeMesh(667, 1001, 15);
        isTriangleSoup = false;
        gdeg = 5;
        edges = [];
        edgeX = [];
        faceX = [];
    end
    
    if nargin < 3
        isTriangleSoup = false;
    end
    mesh.X = X;
    mesh.T = T;
    mesh.nX = size(X,1); nX = mesh.nX;
    mesh.nT = size(T,1); nT = mesh.nT;
    mesh.gdeg = gdeg;
    
    % linear interpolation to form edge and face vertices.
    if sum(size(edgeX))==0
        allEdges = [T(:,[1 2]); T(:,[2 3]); T(:,[3 1])];
        [mesh.edges, ia, ic] = unique(sort(allEdges,2),'rows');
        
        edgeWs = linspace(0,1,gdeg+1);
        intEdgeWs = reshape(edgeWs(2:end-1),1,1,[]);
        mesh.edgeX = mesh.X(mesh.edges(:,1),:).*intEdgeWs + mesh.X(mesh.edges(:,2),:).*(1-intEdgeWs); % nE d nEC
    else
        % when providing edgeX, one must provide edges as well to preserve consistent ordering.
        assert(numel(edges)~=0);
        mesh.edgeX = edgeX;
        mesh.edges = edges;
    end
    mesh.nE = size(mesh.edges,1);
    
    mesh.triangleEdges = permute(reshape(mesh.T(:,[1 2 2 3 3 1]),nT,2,3),[1 3 2]);
    [isF, Find] = ismember(reshape(mesh.triangleEdges,[],2), mesh.edges,'rows');
    [isB, Bind] = ismember(reshape(mesh.triangleEdges,[],2), mesh.edges(:,[2 1]),'rows');
    assert(all(xor(isF, isB)));
    mesh.triangles2edges = reshape(Find + Bind,nT,3);
    mesh.triangles2edges_isForward = reshape(isF,nT,3);
    
    if sum(size(faceX))==0
        [ws, interiorInds] = getBarycentricSamplingWeights(gdeg+1, 1);
        interiorws = reshape(ws(interiorInds,:),1,1,[],3);
        mesh.faceX = mesh.X(mesh.T(:,1),:).*interiorws(1,1,:,1) + mesh.X(mesh.T(:,2),:).*interiorws(1,1,:,2) + mesh.X(mesh.T(:,3),:).*interiorws(1,1,:,3);
    else
        mesh.faceX = faceX;
    end
    mesh.triangleControlPoints = XedgeXfaceXtoControlPoints(mesh);
    mesh.triMaps = poly_triangle_map(permute(mesh.triangleControlPoints,[3 2 1]));
    
    if nargin == 0
        evalVals = mesh.triMaps(getBarycentricSamplingWeights(gdeg+30));
        figure; hold all; axis equal; 
        ti = randi(mesh.nT);
        patch('faces',mesh.T(ti,:),'vertices',mesh.X,'facecolor','green');
        scatter(evalVals(:,ti,1),evalVals(:,ti,2),'k','filled')

        figure; hold all; axis equal; 
        scatter(reshape(mesh.triangleControlPoints(:,1,:),[],1),reshape(mesh.triangleControlPoints(:,2,:),[],1),'c')
        scatter(reshape(mesh.faceX(:,1,:),[],1),reshape(mesh.faceX(:,2,:),[],1),'r.')
        scatter(reshape(mesh.edgeX(:,1,:),[],1),reshape(mesh.edgeX(:,2,:),[],1),'g.')
        scatter(reshape(mesh.X(:,1,:),[],1),reshape(mesh.X(:,2,:),[],1),'b.')
    end
    
%     v1 = X(T(:,1),:);  v2 = X(T(:,2),:);  v3 = X(T(:,3),:);
%     e1 = v2-v1;
%     e2 = v3-v2;
%     e3 = v1-v3;
%     mesh.triangleEdgeNormals = zeros(mesh.nT,3,2);
%     mesh.triangleEdgeNormals(:,1,1) = - e1(:,2);
%     mesh.triangleEdgeNormals(:,1,2) = e1(:,1);
%     mesh.triangleEdgeNormals(:,2,1) = - e2(:,2);
%     mesh.triangleEdgeNormals(:,2,2) = e2(:,1);
%     mesh.triangleEdgeNormals(:,3,1) = - e3(:,2);
%     mesh.triangleEdgeNormals(:,3,2) = e3(:,1);
%     mesh.triangleEdgeNormals = -mesh.triangleEdgeNormals ./ vecnorm(mesh.triangleEdgeNormals,2,3);
%     
%     % you might think the normals should be fliped if triangle is inverted, but empirical testing shows significant instability if you do flip it. ultimately just make sure the mesh isn't inverted...
%     % mesh.triangleEdgeNormals = mesh.triangleEdgeNormals.*((mesh.signedTriAreas>0)-.5)*2;
%     
%     mesh.triangleEdgeLengths = [vecnorm(e1,2,2) vecnorm(e2,2,2) vecnorm(e3,2,2)];
%     
    %% boundary handling. Specific to rectangular boundary mesh!
    isBoundaryEdge = accumarray(mesh.triangles2edges(:),ones(3*size(T,1),1))==1;
    mesh.isBoundaryEdge = isBoundaryEdge;
    edgetangs = X(mesh.edges(:,1),:)-X(mesh.edges(:,2),:);
    edgetangs = edgetangs./vecnorm(edgetangs,2,2);
    isXEdge = abs(edgetangs*[0 1]')<.000001;
    isYEdge = abs(edgetangs*[1 0]')<.000001;
    mesh.isXEdge = isXEdge;
    mesh.isYEdge = isYEdge;
    XvertInds = unique(mesh.edges(isBoundaryEdge & isXEdge,:));
    YvertInds = unique(mesh.edges(isBoundaryEdge & isYEdge,:));
    mesh.isXvert = false(size(X,1),1); mesh.isXvert(XvertInds) = true;
    mesh.isYvert = false(size(X,1),1); mesh.isYvert(YvertInds) = true;
    mesh.isInterior = ~mesh.isXvert & ~mesh.isYvert;
    mesh.isBoundaryTriangle = sum(~mesh.isInterior(mesh.T),2)~=0;
    
    %% compute dA_dt
    TX = permute(reshape(X(T',:),3,[],2),[2 3 1]);
    e123 = TX - circshift(TX,-1,3);
    e312 = circshift(e123,-1,3); % get part of e123 orthogonal to e312
    altitudes = e123 - sum(e123.*e312,2)./vecnorm(e312,2,2).^2.*e312;
    altitudes = altitudes./vecnorm(altitudes,2,2);
    mesh.dAdt = reshape(vecnorm(e312,2,2).*altitudes/2,[],6); %nT, 6=(2(xy) x 3(verts)).
    
    %{
    tind = 22;
    vs = X(T(tind,:),:);
    dirs = reshape(mesh.dAdt(tind,:),2,3)';
    figure; hold all; axis equal; 
    quiver(vs(:,1),vs(:,2),dirs(:,1),dirs(:,2))
    patch('vertices',vs,'faces',[1 2 3],'facecolor','green')
    %}
    
    if ~isTriangleSoup
        %% build edges2triangles matrix
        T2Emat = sparse(repmat([1:mesh.nT]',3,1),mesh.triangles2edges(:),ones(3*mesh.nT,1),mesh.nT,mesh.nE);
        boundaryEdges2Tri = T2Emat(:, mesh.isBoundaryEdge)';
        interiorEdges2Tri = T2Emat(:, ~mesh.isBoundaryEdge)';
        [II,JJ] = find(interiorEdges2Tri);
        [~,perm] = sort(II);
        interiorEdges2TriInds = reshape(JJ(perm),2,[])';
        [II, JJ] = find(boundaryEdges2Tri);
        [~,perm] = sort(II);
        boundaryEdges2TriInds = JJ(perm);
        mesh.edges2triangles = zeros(mesh.nE,2);
        mesh.edges2triangles(mesh.isBoundaryEdge,:) = [boundaryEdges2TriInds boundaryEdges2TriInds];
        mesh.edges2triangles(~mesh.isBoundaryEdge,:) = interiorEdges2TriInds;
    end
end
