function [X,T] = initialMesh(img, initialHorizontalSampling)

    width = size(img,2);
    height = size(img,1);

    xvals = linspace(0,width,initialHorizontalSampling);
    eps = .00001;
    xvals(1) = eps;
    xvals(end) = width - eps;
    dx = xvals(2);
    initialVerticalSampling = ceil(height/dx);
    yvals = linspace(0,height,initialVerticalSampling);
    yvals(1) = eps; yvals(end) = height - eps;
    [Xgrid,Ygrid] = meshgrid(xvals,yvals);
    inds = reshape(1:numel(Ygrid),size(Ygrid,1),size(Ygrid,2));
    TL = inds(1:end-1,1:end-1);
    TR = inds(1:end-1,2:end);
    BL = inds(2:end,1:end-1);
    BR = inds(2:end,2:end);
    % clockwise oriented triangles
    T = [TL(:) TR(:) BL(:); BL(:) TR(:) BR(:)];
    X = [Xgrid(:) Ygrid(:)];

end