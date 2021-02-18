
function [queryEval] = evalBasis(bc_query, deg, coeffs, monomials)
    if deg == 0
        queryEval = ones(size(bc_query,1),1);
        return;
    end
    
    if ~exist('coeffs','var')
        [coeffs, monomials] = getBasisCoeffs(deg);
    end; assert(exist('monomials','var')==1);
    
    u = bc_query(:,1);
    v = bc_query(:,2);
    w = bc_query(:,3);
    submon = double(subs(monomials));
    if numel(submon)==1
        % degree 0 case. monomials doesn't rely on u,v,w so evaluation size is wrong.
        submon = ones(size(bc_query,1),1);
    end
    queryEval = double(submon)*coeffs;
end