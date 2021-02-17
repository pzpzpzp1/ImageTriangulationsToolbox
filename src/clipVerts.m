function X = clipVerts(X, width, height)
    size0 = size(X);
    X = reshape(X,[],2);
    eps = .00001; 
    X(X<eps) = eps;
    X(X(:,1) > width-eps, 1) = width-eps;
    X(X(:,2) > height-eps, 2) = height-eps;
    X = reshape(X,size0);
end