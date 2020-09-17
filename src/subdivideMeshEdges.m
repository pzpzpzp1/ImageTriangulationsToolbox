function [Xnew, Tnew, dividedEdgeInds] = subdivideMeshEdges(mesh, edgeInds, img, n1D)
    %% part 1: find subset of edgeInds that don't conflict with each other. prioritize edges at top of list.
    isCoveredTriangles = false(mesh.nT,1);
    keepEdgeInd = false(numel(edgeInds),1);
    for i=1:numel(edgeInds)
        eind = edgeInds(i);
        ts = mesh.edges2triangles(eind,:);
        if ~any(isCoveredTriangles(ts))
            isCoveredTriangles(ts) = true;
            keepEdgeInd(i) = true;
        end
    end
    edgeInds = edgeInds(keepEdgeInd);
    dividedEdgeInds = edgeInds;
    ne = numel(edgeInds);
    
    %% part 1.5: figure out on each edge where to split it. 50/50 might not be ideal.
    nE = mesh.nE; nT = mesh.nT; nX = mesh.nX;
    T = mesh.T; X = mesh.X; ne = numel(edgeInds);
    if nargin >= 3
        % This doesn't actually seem that effective. Oh well, it can't hurt.
        edgeWs = linspace(1/4,3/4,n1D);
%         edgeWs = linspace(0,1,n1D);
        
        dw = (edgeWs(2)-edgeWs(1))/2;
        X1 = reshape(mesh.X(mesh.edges(edgeInds,1),:),ne,1,2);
        X2 = reshape(mesh.X(mesh.edges(edgeInds,2),:),ne,1,2);
        Xsamples = X1.*(1-edgeWs) + X2.*edgeWs;
        fsamples = double(sampleImage(img,Xsamples)); % [ne n1D 3]
        jumps = vecnorm(fsamples(:,2:n1D,:) - fsamples(:,1:n1D-1,:),2,3); % [ne, n1D-1]
        [~,inds] = max(jumps,[],2);
        ws = 1-(edgeWs(inds) + dw)';
    else
        ws = .5*ones(numel(edgeInds),1);
    end
    
    %% part 2: create subdivided mesh.
    v123 = T(mesh.edges2triangles(edgeInds,1),:);
    v234 = T(mesh.edges2triangles(edgeInds,2),:);
    v23 = mesh.edges(edgeInds,:);
    [II,JJ] = find(~squeeze(sum(reshape(v123,ne ,1,3) == v23,2)));
    [~,perm] = sort(II);
    v1 = v123(sub2ind(size(v123),[1:ne ]',JJ(perm)));
    [II,JJ] = find(~squeeze(sum(reshape(v234,ne ,1,3) == v23,2)));
    [~,perm] = sort(II);
    v4 = v234(sub2ind(size(v234),[1:ne ]',JJ(perm)));
    v2 = v23(:,1);
    v3 = v23(:,2);
    v5 = ((nX+1):(nX+ne))';
    Xnew = [X; ws.*X(mesh.edges(edgeInds,1),:) + (1-ws).*X(mesh.edges(edgeInds,2),:)];
    Tnew = [T(~isCoveredTriangles,:); v1 v2 v5; v5 v2 v4; v3 v5 v4; v1 v5 v3];
    Tnew = unique(sort(Tnew,2),'rows'); % get rid of any duplicated triangles from splitting boundary edge.
    
    % reorient inverted triangles
    v1 = Xnew(Tnew(:,1),:);  v2 = Xnew(Tnew(:,2),:);  v3 = Xnew(Tnew(:,3),:);
    e12 = [v1-v2]; e12(1,3) = 0; e23 = [v2-v3]; e23(1,3) = 0;
    flippedTris = (cross(e12,e23)/2)*[0 0 1]' < 0;
    Tnew(flippedTris,:) = Tnew(flippedTris,[1 3 2]);
    
    %{
    figure; 
    subplot(1,2,2); image(img); hold all; axis equal; axis off;
    patch('vertices',mesh.X,'faces',mesh.edges(dividedEdgeInds,[1 2 1]),'facecolor','none','linewidth',2,'edgecolor','c');
    patch('vertices',Xnew,'faces',Tnew,'facecolor','green','facealpha',.1);
    subplot(1,2,1); image(img); hold all; axis equal; axis off;
    patch('vertices',mesh.X,'faces',mesh.T,'facecolor','green','facealpha',0,'edgecolor','c');
    %}
end