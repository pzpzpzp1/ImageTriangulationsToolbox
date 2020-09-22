function X = clipVerts(X, width, height)
    eps = .00001; 
    X(X<eps) = eps;
    X(X(:,1) > width-eps, 1) = width-eps;
    X(X(:,2) > height-eps, 2) = height-eps;
end