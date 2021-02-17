function [w, grad] = RegularSimplexToBarycentric(u)
    if nargin==0
        N=100;
        u = abs(rand(N,2)); 
        u(u(:,1)+u(:,2)>1,:)=[];
        % scatter(u(:,1),u(:,2))
    end
    N = size(u,1);
    p = [[0 0; 1 0; 0 1]'; 1 1 1];
    u1 = u; u1(:,3)=1;
    invp = inv(p);
    w = (invp*u1')';
    grad = invp(:,1:2);
    
    if nargin==0
        % verify fdiff
        eps = 1e-6;
        pert = randn(size(u));
        wP = RegularSimplexToBarycentric(u+eps*pert);
        wM = RegularSimplexToBarycentric(u-eps*pert);
        fdiff = (wP-wM)/(2*eps);
        adiff = (grad*pert')';
        assert(norm(fdiff-adiff)<eps);
    end
end