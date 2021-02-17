function [u, grad] = barycentricToRegularSimplex(w)
    if nargin==0
        N=100;
        w = abs(randn(N,3)); w=w./sum(w,2);
    end
    N = size(w,1);
    p = [0 0; 1 0; 0 1];
    u = w*p;
    grad = p';
end