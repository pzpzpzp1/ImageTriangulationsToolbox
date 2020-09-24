files = dir('../images/*.jpg');

for i=2:numel(files)
    fname = files(i).name;
    fullname = [files(i).folder '\' fname];
    img = imread(fullname);
    prevsal = [];
    if exist(fname,'file')
        prevsal = imread(fname);
    end
    
    boost = getBoostMap(img, prevsal);
    
    if norm(boost,'fro')~=0
        imwrite(boost, fname);
    end
end