% Computes colors on vertices per triangle in two different ways. 
% First using all sampled colors over the triangle. 
% Second using a subset of sampled colors dicated by interiorInds.
% The second may be better for rendering purposes. Less sensitive to edge artifacts.
% mesh must have X and T
% sampleVals             [nT n 3]
% interiorInds           [n]
function [colors, colorsAlt] = constantGetColorsFromSamples(sampleVals, interiorInds, saliencySamples)
    %% compute full mean per triangle
    normalizer = sum(saliencySamples,1);
    colors = squeeze(sum(sampleVals.*saliencySamples,1)./normalizer);
    
    %% compute full mean per triangle
    normalizerAlt = sum(saliencySamples(interiorInds,:),1);
    colorsAlt = squeeze(sum(sampleVals(interiorInds,:,:).*saliencySamples(interiorInds,:),1)./normalizerAlt);
    
end