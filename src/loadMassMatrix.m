% computes int_T phi_i phi_j where T is a standard triangle:(0,0),(1,0),(0,1).
% and phi_i, phi_j are the FEM basis funcs of a certain degree on that triangle.
function [K, Ki] = loadMassMatrix(degree)
    persistent cachedKs cachedKis;
    if isempty(cachedKs) || degree > numel(cachedKs) || numel(cachedKs{degree})==0
        if degree == 0
            error('you shouldnt need to load the constant case mass matrix');
        elseif degree == 1
            cachedKs{degree} = ones(3)/24; cachedKs{degree}(1:4:end) = 1/12;
            cachedKis{degree} = (4*eye(3) - 1)*6;
        else
            error('unhandled mass matrix');
        end
    end
    K = cachedKs{degree};
    Ki = cachedKis{degree};
end