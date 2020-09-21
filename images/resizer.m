files = dir('./*.jpg');

for i=1:numel(files)
    img = imread(files(i).name);
    if any(size(img,[1 2])>1000)
        factor = max(size(img,[1 2])/1000);
        img2 = imresize(img,1/factor);
        imwrite(img2, [files(i).name]);
    end
end

