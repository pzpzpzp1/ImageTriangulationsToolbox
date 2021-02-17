
function [queryEval] = evalBasis(bc_query, deg, coeffs, monomials)
    if ~exist('coeffs','var')
        [coeffs, monomials] = getBasisCoeffs(deg);
    end; assert(exist('monomials','var')==1);
    
    u = bc_query(:,1);
    v = bc_query(:,2);
    w = bc_query(:,3);
    queryEval = double(subs(monomials))*coeffs;
end