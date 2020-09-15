width = 255;
height = 255;
imdata = zeros(height,width,3);
xvals = 1:width;
yvals = 1:height;
[X,Y] = meshgrid(xvals,yvals);

for i=1:10
    direc = abs(randn(2,1)); direc = direc/norm(direc);
%     direc = [1 0];
    imdata(:,:,1) = 100*rand*(X.*direc(1)+Y.*direc(2))+rand*100;
    imdata(:,:,2) = 100*rand*(X.*direc(1)+Y.*direc(2))+rand*100;
    imdata(:,:,3) = 100*rand*(X.*direc(1)+Y.*direc(2))+rand*100;
    imdata = imdata/max(imdata,[],'all')*256;
    
    imdata(:,:,1) = repmat(1:255,255,1);
    imdata(:,:,2) = repmat(1:255,255,1);
    imdata(:,:,3) = repmat(1:255,255,1);
    
%     K = mean(imdata,3);
%     imdata(:,:,1) = K;
%     imdata(:,:,2) = K;
%     imdata(:,:,3) = K;
    figure;image(uint8(imdata)); axis equal;

    imwrite(uint8(imdata),sprintf('gradient%d.PNG',i));
%     imwrite(uint8(imdata),sprintf('gradientHorizontalGray.PNG',i));
%     imwrite(uint8(imdata),sprintf('gradientVerticalGray.PNG',i));
%     imwrite(uint8(imdata),sprintf('gradientDiag.PNG',i));

end