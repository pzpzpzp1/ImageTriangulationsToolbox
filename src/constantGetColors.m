function colors = constantGetColors(img, X, T, integral1DNsamples)
    ws = getBarycentricSamplingWeights(integral1DNsamples);
    samplePoints = getSamplePointsFromBarycentricWeights(ws, X, T);
    n = size(ws,1);
    nT = size(T,1);
    
    xind = ceil(reshape(samplePoints(:,:,1),[],1));
    yind = ceil(reshape(samplePoints(:,:,2),[],1));
    linearInds = sub2ind(size(img,[1 2]), yind, xind);
    flatimg = reshape(img,[],3);
    sampleVals = reshape(flatimg(linearInds,:), n, nT, 3);
    colors = ceil(squeeze(mean(sampleVals,1)));
end