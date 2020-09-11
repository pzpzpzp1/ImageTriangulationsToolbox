% size of output:
% 1 placeholder for number of samples per edge, 
% 3 edges per triangle, 
% 6 perturbbations dofs for 3 vertices, 
% 2 xy coords for output response.
function [zf, zb] = loadVSamplingMatrix
    persistent cachedZf cachedZb;
    if isempty(cachedZf) || isempty(cachedZb) 
        cachedZf = zeros(1, 3, 6, 2);
        cachedZf(1, [1], 1, 1) = 1; % vertex 1 is part of edge 3 and 1.
        cachedZf(1, [1], 2, 2) = 1; % vertex 1 is part of edge 3 and 1.
        cachedZf(1, [2], 3, 1) = 1; % vertex 2 is part of edge 1 and 2.
        cachedZf(1, [2], 4, 2) = 1; % vertex 2 is part of edge 1 and 2.
        cachedZf(1, [3], 5, 1) = 1; % vertex 3 is part of edge 2 and 3.
        cachedZf(1, [3], 6, 2) = 1; % vertex 3 is part of edge 2 and 3.
        
        cachedZb = zeros(1, 3, 6, 2);
        cachedZb(1, [3], 1, 1) = 1; % vertex 1 is part of edge 3 and 1.
        cachedZb(1, [3], 2, 2) = 1; % vertex 1 is part of edge 3 and 1.
        cachedZb(1, [1], 3, 1) = 1; % vertex 2 is part of edge 1 and 2.
        cachedZb(1, [1], 4, 2) = 1; % vertex 2 is part of edge 1 and 2.
        cachedZb(1, [2], 5, 1) = 1; % vertex 3 is part of edge 2 and 3.
        cachedZb(1, [2], 6, 2) = 1; % vertex 3 is part of edge 2 and 3.
    end
    zf = cachedZf;
    zb = cachedZb;
end