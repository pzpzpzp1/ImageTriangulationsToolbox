% size of output:
% 1 placeholder for number of samples per edge, 
% 3 edges per triangle, 
% 6 perturbbations dofs for 3 vertices, 
% 2 xy coords for output response.
function z = loadVSamplingMatrix
    persistent cachedZ;
    if isempty(cachedZ) 
        cachedZ = zeros(1, 3, 6, 2);
        cachedZ(1, [3 1], 1, 1) = 1; % vertex 1 is part of edge 3 and 1.
        cachedZ(1, [3 1], 2, 2) = 1; % vertex 1 is part of edge 3 and 1.
        cachedZ(1, [1 2], 3, 1) = 1; % vertex 2 is part of edge 1 and 2.
        cachedZ(1, [1 2], 4, 2) = 1; % vertex 2 is part of edge 1 and 2.
        cachedZ(1, [2 3], 5, 1) = 1; % vertex 3 is part of edge 2 and 3.
        cachedZ(1, [2 3], 6, 2) = 1; % vertex 3 is part of edge 2 and 3.
    end
    z = cachedZ;
end