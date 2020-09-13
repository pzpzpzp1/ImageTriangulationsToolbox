function [X,T] = initialHexLatticeMesh(width, height, initialHorizontalSampling)
    assert(initialHorizontalSampling >= 3);
    
    xvals = linspace(0,width,initialHorizontalSampling);
    eps = .00001;
    xvals(1) = eps;
    xvals(end) = width - eps;
    dx = xvals(2);
    
    initialVerticalSampling = ceil(height/dx/sin(pi/3));
    yvals = linspace(0,height,initialVerticalSampling);
    yvals(1) = eps; yvals(end) = height - eps;
    [Xgrid, Ygrid] = meshgrid(xvals,fliplr(yvals));
    Xgrid(2:2:end,:) = Xgrid(2:2:end,:) + dx/2;
    
    inds = reshape(1:numel(Ygrid),size(Ygrid,1),size(Ygrid,2));
    TL = inds(1:2:end-1,1:end-1);
    TR = inds(1:2:end-1,2:end);
    BL = inds(2:2:end,1:end-1);
    BR = inds(2:2:end,2:end);
    T1 = [TL(:) TR(:) BL(:); BL(:) TR(:) BR(:)];
    TL2 = inds(2:2:end-1,1:end-1);
    TR2 = inds(2:2:end-1,2:end);
    BL2 = inds(3:2:end,1:end-1);
    BR2 = inds(3:2:end,2:end);
    T2 = [TL2(:) TR2(:) BR2(:); BL2(:) TL2(:) BR2(:)];
    Tm = [T1;T2];
    
    padX = Xgrid(2:2:end, 1)*0+eps;
    padY = Ygrid(2:2:end, 1);
    padInds = (1:numel(padX)) + numel(Xgrid);
    padUind = inds(1:2:end,1);
    padDind = inds(3:2:end,1);
    padRind = inds(2:2:end,1);
    padT1 = [padInds(:) padUind(1:numel(padInds)) padRind(1:numel(padInds));];
    padT2 = [padInds(1:numel(padDind))' padRind(1:numel(padDind)) padDind(:)];
    T = [Tm; padT1; padT2];
    T = T(:,[1 3 2]);
    X = [Xgrid(:) Ygrid(:); padX(:) padY(:)];
    X(X(:,1) > width) = width - eps;
    
    %{
    figure; hold all; axis equal; scatter(Xgrid(:), Ygrid(:),'r.')
    patch('faces',Tm,'vertices',X,'facecolor','green','facealpha',.5);
    patch('faces',padT1,'vertices',X,'facecolor','blue','facealpha',.5);
    patch('faces',padT2,'vertices',X,'facecolor','red','facealpha',.5);
    %}
end