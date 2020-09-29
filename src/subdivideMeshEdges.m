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
    [II,JJ] = find(~permute(sum(reshape(v123,ne ,1,3) == v23,2),[1 3 2]));
    [~,perm] = sort(II);
    v1 = v123(sub2ind(size(v123),[1:ne ]',JJ(perm)));
    ip1 = mod(JJ(perm)+1,3); ip1(ip1==0)=3;
    ip2 = mod(JJ(perm)+2,3); ip2(ip2==0)=3;
    v2 = v123(sub2ind(size(v123),[1:ne ]',ip1));
    v3 = v123(sub2ind(size(v123),[1:ne ]',ip2));
    [II,JJ] = find(~permute(sum(reshape(v234,ne ,1,3) == v23,2),[1 3 2]));
    [~,perm] = sort(II);
    v4 = v234(sub2ind(size(v234),[1:ne ]',JJ(perm)));
    v5 = ((nX+1):(nX+ne))';
    Xnew = [X; ws.*X(mesh.edges(edgeInds,1),:) + (1-ws).*X(mesh.edges(edgeInds,2),:)];
    dedupBoundary = ~mesh.isBoundaryEdge(edgeInds);
    Tnew = [T(~isCoveredTriangles,:); v1 v2 v5; v1 v5 v3; v5(dedupBoundary) v2(dedupBoundary) v4(dedupBoundary); v3(dedupBoundary) v5(dedupBoundary) v4(dedupBoundary); ];
    assert(size(unique(sort(Tnew,2),'rows'),1) == size(Tnew,1)); % should be no dupes
    
    %{
    figure; 
    subplot(1,2,2); image(img); hold all; axis equal; axis off;
    patch('vertices',mesh.X,'faces',mesh.edges(dividedEdgeInds,[1 2 1]),'facecolor','none','linewidth',2,'edgecolor','c');
    patch('vertices',Xnew,'faces',Tnew,'facecolor','green','facealpha',.1);
    subplot(1,2,1); image(img); hold all; axis equal; axis off;
    patch('vertices',mesh.X,'faces',mesh.T,'facecolor','green','facealpha',0,'edgecolor','c');

    figure; 
    ax1=subplot(1,2,2); hold all; axis equal; axis off;
    patch('vertices',mesh.X,'faces',mesh.edges(dividedEdgeInds,[1 2 1]),'facecolor','none','linewidth',2,'edgecolor','c');
    patch('vertices',Xnew,'faces',Tnew,'facecolor','blue','facealpha',.1);
    patch('vertices',mesh.X,'faces',mesh.T(unique(mesh.edges2triangles(dividedEdgeInds,:)),:),'facecolor','red','facealpha',.5);
    ax2=subplot(1,2,1); hold all; axis equal; axis off;
    patch('vertices',mesh.X,'faces',mesh.T,'facecolor','blue','facealpha',.1);
    patch('vertices',mesh.X,'faces',mesh.edges(dividedEdgeInds,[1 2 1]),'facecolor','none','linewidth',2,'edgecolor','c');  
    ax1.Clipping='off';ax2.Clipping='off';
    %}
end