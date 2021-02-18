function [out, outmonomials, outdiffmonomials] = getBasisCoeffs(n)
    assert(n >= 0);
    if n==0
        out = 1;
        outmonomials = 1;
        outdiffmonomials = [0 0 0];
        return;
    end
    persistent cachedWs cachedMons cachedDiffMons
    if isempty(cachedWs) || n > numel(cachedWs) || numel(cachedWs{n})==0
        ws = getBarycentricSamplingWeights(n+1); nC = size(ws,1); % barycentric control points
        syms u v w real; sympoly = expand((u+v+w)^n);
        [~,monomials] = coeffs(sympoly); nm = numel(monomials);
        diffmon = [diff(monomials',u) diff(monomials',v) diff(monomials',w)];
        u = ws(:,1)'; v = ws(:,2)'; w = ws(:,3)';
        outcoeffs = inv(double(subs(monomials'))');
        
        cachedWs{n} = outcoeffs;
        cachedMons{n} = monomials;
        cachedDiffMons{n} = diffmon;
    end    
    out = cachedWs{n};
    outmonomials = cachedMons{n};
    outdiffmonomials = cachedDiffMons{n};
end
