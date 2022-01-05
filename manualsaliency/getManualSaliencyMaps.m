files = dir('../images/*.jpg');
whitelist = {'jinx','vi','silco','jinx2','jinx3'};
whitelist = {'blu','goldenfrog','frog','graha'};
whitelist = {'frog'};

for i=1:numel(files)
    fname = files(i).name;
    [~,name,ext] = fileparts(files(i).name);
    if exist('whitelist','var') && ~any(strcmp(whitelist,name))
        continue;
    end
    
    fullname = [files(i).folder '\' fname];
    img = imread(fullname);
    prevsal = [];
    if exist(['./' fname],'file')
        prevsal = imread(fname);
    end
    
    boost = getBoostMap(img, prevsal);
    
    if norm(boost,'fro')~=0
        imwrite(boost, fname);
    end
end