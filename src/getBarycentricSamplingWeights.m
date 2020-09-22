
function [ws, interiorInds] = getBarycentricSamplingWeights(n)
    assert(n >= 5);
    persistent cachedWs cachedInteriorInds
    if isempty(cachedWs) || n > numel(cachedWs) || numel(cachedWs{n})==0
        % n would need to be 1e6 for this to become floating point imprecise.
        eps = 1e-6;
        
        u = linspace(0,1,n);
        v = linspace(0,1,n);
        [U,V] = ndgrid(u,v);
        W = 1-U-V;
        
        keep = find(W>=-eps);
        assert(numel(keep)==n*(n+1)/2); % another thresholding sanity check
        cachedWs{n} = [U(keep) V(keep) W(keep)];
        layersdeep = 3;
        cachedInteriorInds{n} = ~any(abs(cachedWs{n}) < (layersdeep-1)/(n-1) + eps,2); % 1 layer interior
    end
    ws = cachedWs{n};
    interiorInds = cachedInteriorInds{n};
end

 %% sanity check visualization
    %{
    X = [0 0;1 0;0 1;]; T=[1 2 3];
    samplePoints1 = ws(interiorInds,:)*X;
    samplePoints2 = ws(~interiorInds,:)*X;
    figure; hold all; axis equal;
    scatter(samplePoints2(:,1),samplePoints2(:,2),20,'r','filled')
    scatter(samplePoints1(:,1),samplePoints1(:,2),30,'b')
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