
function [ws, interiorInds, tris] = getBarycentricSamplingWeights(n, layersdeep)
    if nargin == 0
        clear all; close all;
        n=15;
        layersdeep = 3;
    end
    
    assert(n >= 2);
    persistent cachedWs cachedInteriorInds cachedTris
    if isempty(cachedWs) || n > numel(cachedWs) || numel(cachedWs{n})==0
        % n would need to be 1e6 for this to become floating point imprecise.
        eps = 1e-6;
        
        u = linspace(0,1,n);
        v = linspace(0,1,n);
        [U,V] = ndgrid(u,v);
        W = 1-U-V;
        
        keep = find(W>=-eps);
        assert(numel(keep)==n*(n+1)/2); % another thresholding sanity check
        
        [cachedWs{n}, perm] = sortrows([U(keep) V(keep) W(keep)]);
        cachedWs{n}(abs(cachedWs{n}) < eps) = 0;
        if ~exist('layersdeep','var'); layersdeep = 2; end;
        cachedInteriorInds{n} = ~any(abs(cachedWs{n}) < (layersdeep-1)/(n-1) + eps,2); % 1 layer interior
        cachedInteriorInds{n}(abs(cachedInteriorInds{n}) < eps) = 0;
        
        % get triangulation
        notkeep = W < -eps;
        Linds = zeros(size(W)); Linds(:)=1:numel(Linds);
        t1 = Linds(1:end-1,1:end-1); t2 = Linds(2:end,1:end-1); t3 = Linds(2:end,2:end); t4 = Linds(1:end-1,2:end);
        fulltris = [t1(:) t2(:) t4(:); t4(:) t2(:) t3(:)];
        fulltris(any(ismember(fulltris,Linds(notkeep)),2),:)=[]; % remove half of the triangles
        inds = knnsearch(cachedWs{n},[U(:) V(:) W(:)]); % lazy re-indexing
        fulltris = inds(fulltris);
        cachedTris{n} = fulltris;
    end
    ws = cachedWs{n};
    interiorInds = cachedInteriorInds{n};
    tris = cachedTris{n};
end

 %% sanity check visualization
    %{
    onelayerin = ~any(abs(cachedWs{n}) < (1-1)/(n-1) + eps,2); % 1 layer interior
    X = [0 0;1 0;0 1;]; T=[1 2 3];
    samplePoints1 = ws(interiorInds,:)*X;
    samplePoints2 = ws(~interiorInds,:)*X;
    samplePoints3 = ws(~onelayerin,:)*X;

    figure; hold all; axis equal;
    plot([0 1 0 0],[0 0 1 0],'k','linewidth',2);
    scatter(samplePoints2(:,1),samplePoints2(:,2),30,'r','filled')
    scatter(samplePoints1(:,1),samplePoints1(:,2),30,'r','filled')
    scatter(samplePoints3(:,1),samplePoints3(:,2),30,'g','filled')
 
    axis off; set(gcf,'color','w'); 
    
    %}
 
    %{
    [K, Ki] = loadMassMatrix(1);
    X = [0 0;1 0;0 1;]; T=[1 2 3];
    samplePoints = ws*X;
    figure; 
    subplot(1,3,1);hold all; axis equal;
    scatter(samplePoints(:,1),samplePoints(:,2),20,ws(:,1),'filled')
    colorbar;
    subplot(1,3,2);hold all; axis equal;
    scatter(samplePoints(:,1),samplePoints(:,2),20,ws(:,2),'filled')
    colorbar;
    subplot(1,3,3);hold all; axis equal;
    scatter(samplePoints(:,1),samplePoints(:,2),20,ws(:,3),'filled')
    colorbar;
    Ksamp = ws'*ws/(2*size(ws,1));
    norm(Ksamp - K)
    %}  