% phi_edges_oneTri: (samples per edge) (3 edges) (3 phis)
function z = getEdgeBarycentricSamplingWeights(n1D)
    assert(n1D >= 3);
    persistent cachedWs;
    if isempty(cachedWs) || n1D > numel(cachedWs) || numel(cachedWs{n1D})==0
        
        % triangle: [v1 e1 v2 e2 v3 e3 goto_start]
        edgeWs = linspace(0,1,n1D);
        phi_edges_oneTri = zeros(n1D,3,3);
        phi_edges_oneTri(:,1,1) = fliplr(edgeWs);
        phi_edges_oneTri(:,2,2) = fliplr(edgeWs);
        phi_edges_oneTri(:,3,3) = fliplr(edgeWs);
        phi_edges_oneTri(:,3,1) = edgeWs;
        phi_edges_oneTri(:,1,2) = edgeWs;
        phi_edges_oneTri(:,2,3) = edgeWs;
        cachedWs{n1D} = phi_edges_oneTri;
    end
    z = cachedWs{n1D};
end