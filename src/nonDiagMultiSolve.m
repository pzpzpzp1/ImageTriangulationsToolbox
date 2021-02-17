% solves many small linear systems in one go. 
% As:n x n2 x m 
% bs: n x 1 x m
% solutions: n2 x m
%%
function solutions = nonDiagMultiSolve(As,bs)

if nargin == 0
    n = 465; n2=21; m=3000; K = 3;
%     n = 120; n2=6; m=3000; K = 3;
    n = 1; n2=3; m=54; K = 3;
    As = randn(n,n2,m);
    bs = randn(n,K,m);
end

n = size(As,1);
n2 = size(As,2);
m = size(As,3);
assert(size(bs,1)==n);
K = size(bs,2);
assert(size(bs,3)==m);


% if it's not fairly sparse perfectly diagonal blocks, the sparse block diag version is slower than just for looping.
% tic
% solutions = zeros(n2,K,m);
% for i = 1:m
%     for j=1:K
%         solutions(:,j,i) = As(:,:,i)\bs(:,j,i);
%     end
% end
% toc

% gpu batchop is at worst 2x slower than for loop, and at best 10x faster for the expected range of input dimensions. 
% for simplicity, lets go with batchop for all cases.
% When K=3, batchop is 5x faster on full range of expected input dimensions.
% tic
try
    solutions = gather(batchop('leastsq', gpuArray(As), gpuArray(bs)));
catch
    % n=1 case doesn't work on batchop for some reason.
    solutions = zeros(n2,K,m);
    for i = 1:m
        for j=1:K
            solutions(:,j,i) = As(:,:,i)\bs(:,j,i);
        end
    end
end
% toc 

% norm(solutions(:)-solutions2(:))


end