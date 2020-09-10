function ws = getBarycentricSamplingWeights(n)
    assert(n >= 3 && n <= 100) % some reasonable bounds
    
    persistent cachedWs;
    if isempty(cachedWs) || n > numel(cachedWs) || numel(cachedWs{n})==0
        u = linspace(0,1,n);
        v = linspace(0,1,n);
        [U,V] = ndgrid(u,v);
        W = 1-U-V;
        keep = find(W>=0);
        cachedWs{n} = [U(keep) V(keep) W(keep)];
    end
    ws = cachedWs{n};
end