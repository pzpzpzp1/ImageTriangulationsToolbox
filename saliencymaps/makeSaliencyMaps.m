clear all; close all;
files = dir('./../images/*.jpg');

figure;
params = defaultSaliencyParams;
for i=1:numel(files)
    img = initializeImage([files(i).folder '/' files(i).name]);
    
    salmap = makeSaliencyMap(img,params);
    bigMap = imresize(salmap.data,img.size(1:2));
    
    imwrite(bigMap, files(i).name);
    
    
    fname = [files(i).folder '/' files(i).name];
    im = imread(fname);
    clf; image(repmat(rgb2gray(im),1,1,3)); hold all; axis equal;
    imh = imagesc(bigMap); hold all; axis equal;
    imh.AlphaData = .5;
    pause;
end


