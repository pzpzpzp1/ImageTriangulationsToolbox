% Computes colors on vertices per triangle in two different ways. 
% First using all sampled colors over the triangle. 
% Second using a subset of sampled colors dicated by interiorInds.
% The second may be better for rendering purposes. Less sensitive to edge artifacts.
% mesh must have X and T
% sampleVals             [nT n 3]
% interiorInds           [n]
function [colors, colorsAlt] = constantGetColorsFromSamples(sampleVals, interiorInds)
    %% compute full mean per triangle
    colors = squeeze(mean(sampleVals,1));
    
    %% compute full mean per triangle
    colorsAlt = squeeze(mean(sampleVals(interiorInds,:,:),1));

end