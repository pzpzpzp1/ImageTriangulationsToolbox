function [X,T] = initialGridMesh(width, height, initialHorizontalSampling, randslant)
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
    slantDirection = rand(size(BR))>.5;
    % clockwise oriented triangles
    if randslant
        T1 = [TL(slantDirection) TR(slantDirection) BL(slantDirection); BL(slantDirection) TR(slantDirection) BR(slantDirection)];
        T2 = [TL(~slantDirection) BR(~slantDirection) BL(~slantDirection); BR(~slantDirection) TL(~slantDirection) TR(~slantDirection)];
        T = [T1;T2];
    else
        T = [TL(:) TR(:) BL(:); BL(:) TR(:) BR(:)];
    end
    X = [Xgrid(:) Ygrid(:)];

end