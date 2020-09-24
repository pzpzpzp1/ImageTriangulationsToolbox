clear all; close all;
files = dir('../images/*.jpg');

qualities = [0];
for i=1:numel(files)
    srcname = [files(i).folder '/' files(i).name];
    img = imread(srcname); 
    [~,name,ext]=fileparts(files(i).name);
    
    for j=1:numel(qualities)
        quality = qualities(j);
        outname = [name '_' num2str(quality) '.jpg'];
        imwrite(img,outname,'Quality',quality);
    end
end

