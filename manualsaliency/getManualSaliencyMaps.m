files = dir('../images/*.jpg');

for i=2:numel(files)
    fname = files(i).name;
    fullname = [files(i).folder '\' fname];
    img = imread(fullname);
    boost = getBoostMap(img);
    imwrite(boost, fname);
end