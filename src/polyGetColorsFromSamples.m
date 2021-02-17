% Computes colors on control points per triangle in two different ways. 
% First using all sampled colors over the triangle. 
% Second using a subset of sampled colors dictated by interiorInds.
% The second may be better for rendering purposes. Less sensitive to edge artifacts.
function colorCPs = polyGetColorsFromSamples(ws, f_triangle, sal_triangle, interiorInds, polyparams, geodets, colstrat)
    if ~exist('colstrat','var')
        colstrat = colorStrategy.full;
    end
    nD = size(f_triangle,1); % num data per tri
    basisEvals = evalBasis(ws, polyparams.cdeg);
    nB = size(basisEvals,2); % num basis funcs
    metric = sqrt(abs(geodets.*sal_triangle));
    LHS = permute(metric .* reshape(basisEvals,nD,1,nB),[1 3 2]);
    RHS = permute(metric.*f_triangle,[1 3 2]);
    
    if colstrat == colorStrategy.both || colstrat == colorStrategy.int
        colorCPs.int = nonDiagMultiSolve(LHS(interiorInds,:,:), RHS(interiorInds,:,:)); % nB,3,nT
    end
    if colstrat == colorStrategy.both || colstrat == colorStrategy.full
        colorCPs.full = nonDiagMultiSolve(LHS, RHS); 
    end
end