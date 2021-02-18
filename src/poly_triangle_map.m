function [triMaps] = poly_triangle_map(V)
	if nargin == 0
        T = randi(100);
        d=3;
        V = randn(6,d,T);
    end
    % get degree based on size of V. if V doesn't have enough points, error.
    ctrpts = size(V, 1);
	n = (-3 + sqrt(1 + ctrpts * 8))/2;
    assert(round(n)==n);
    
    [coeffs, monomials, diffmon] = getBasisCoeffs(n);
    
    triMaps = @(query) evalTri(query, coeffs, monomials, diffmon, V);
    
    if nargin==0
        [query] = getBarycentricSamplingWeights(10); 
        [~, grad] = triMaps(query);
        
        eps = 1e-6; pert = randn(size(query));
        resultP = triMaps(query + eps*pert);
        resultM = triMaps(query - eps*pert);
        fdiff = (resultP-resultM)/(2*eps);
        adiff = sum(grad.*reshape(pert,[],1,1,3),4);
        assert(norm(fdiff(:) - adiff(:),'fro') < eps);
    end
    
end


function [queryEval, geometricGrad, geodets] = evalTri(bc_query, coeffs, monomials, diffmon, V)
    d = size(V,2); nC = size(V,3); nM = size(monomials,2);
    nQ = size(bc_query,1); u=bc_query(:,1); v=bc_query(:,2); w=bc_query(:,3);
    evaluatedBasisFuncs = evalBasis(bc_query, -1, coeffs, monomials);
    
    % nQ nM 3
    evaluatedDerivBasisFuncs = permute(sum(double(subs(reshape(diffmon,1,nM,3))).*reshape(coeffs,1,nM,1,nM),2),[1 4 3 2]);
    if nM==3 % in the linear case diffmon has no dependence on u,v,w. meaning we need to repmat to make the dimensions right.
        evaluatedDerivBasisFuncs = repmat(evaluatedDerivBasisFuncs,nQ,1,1);
    end
    
    queryEval = reshape(evaluatedBasisFuncs*reshape(V,[],d*nC),nQ,d,nC);
    queryEval = permute(queryEval,[1 3 2]); % [nQ nC d]
    
    if d==2
        % [nQ nM 3]                     [nM d nC] -> [nQ nC d 3]
        diffQueryEval = permute(reshape(sum(evaluatedDerivBasisFuncs .* reshape(V,1,nM,1,d,nC),2),nQ,3,d,nC),[1 4 3 2]);
        
        %[U, UGrad] = barycentricToRegularSimplex(bc_query);
        %[W, Wgrad] = RegularSimplexToBarycentric(U);
        
        Wgrad = [-1 -1; 1 0; 0 1];
        geometricGrad = reshape(sum(diffQueryEval .* reshape(Wgrad,1,1,1,3,d),4),nQ,nC,d,d);

        % get jacobian determinant
        geodets = reshape(get2DDeterminants(reshape(geometricGrad,[],2,2)),nQ,nC); % nQ, nC
    end
end















