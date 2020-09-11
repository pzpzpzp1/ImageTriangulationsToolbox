% takes in image and sample locations. returns sampled values.
% samplePoints must be size (n,m,...,2)
% sampleVals will be size (n,m,...,3)
function sampleVals = sampleImage(img, samplePoints)
    nd = ndims(samplePoints);
    assert(size(samplePoints,nd)==2);
    
    mid = numel(samplePoints)/2;
    xind = ceil(samplePoints(1:mid)');
    yind = ceil(samplePoints((mid+1):end)');
    linearInds = sub2ind(size(img,[1 2]), yind, xind);
    flatimg = reshape(img,[],3);
    sampleVals = zeros([size(samplePoints,1:nd-1),3],'uint8');
    % If there's an error here, then the samplePoint is trying to access locations off of the image. No bounds check for now. Just don't read off the image please.
    sampleVals(:) = flatimg(linearInds,:);
end