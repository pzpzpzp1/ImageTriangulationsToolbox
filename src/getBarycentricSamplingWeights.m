function ws = getBarycentricSamplingWeights(n)
    assert(n >= 3 && n <= 5000) % some reasonable bounds
    
    persistent cachedWs;
    if isempty(cachedWs) || n > numel(cachedWs) || numel(cachedWs{n})==0
        u = linspace(0,1,n);
        v = linspace(0,1,n);
        [U,V] = ndgrid(u,v);
        W = 1-U-V;
        % n would need to be 1e6 for this to become floating point imprecise.
        keep = find(W>=-1e-6);
        assert(numel(keep)==n*(n+1)/2); % another thresholding sanity check
        cachedWs{n} = [U(keep) V(keep) W(keep)];
        
    end
    ws = cachedWs{n};
    
    %% sanity check
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
    
end