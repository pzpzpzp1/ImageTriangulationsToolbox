function [extra, energy, colormaps] = polyComputeEnergyCore(img, mesh, integral1DNsamples, salmap, polyparams, colstrat, colormaps)
    if numel(salmap)==0
        salmap = img(:,:,1)*0+1;
    end
    
    extra = {};
    X = mesh.X; T = mesh.T; nT = size(T,1); 
    % generate sample locations in barycentric coords
    [ws, interiorInds, bcTris] = getBarycentricSamplingWeights(integral1DNsamples); n = size(ws,1); 
    [samplePoints, ~, geodets] = mesh.triMaps(ws); % world space points, and jacobian det
    
    % perform sampling
    f_triangle = double(sampleImage(img, samplePoints)); % image colors
    sal_triangle = double(sampleImage(salmap, samplePoints)); % saliency map values
    
    if ~exist('colormaps','var')
        % get colors per triangle
        colorCPs = polyGetColorsFromSamples(ws, f_triangle, sal_triangle, interiorInds, polyparams, geodets, colstrat);
        extra.colorCPs = colorCPs;
        extra.samplePoints = samplePoints;
        
        % build colormap
        if colstrat == colorStrategy.full || colstrat == colorStrategy.both
            colormaps = poly_triangle_map(colorCPs.full);
        else
            colormaps = poly_triangle_map(colorCPs.int);
        end
    end
    evalColors = colormaps(ws);
    
    % compute energy integrated over triangles
    extra.perTriangleRGBError = reshape(sum((evalColors - f_triangle).^2 .* abs(geodets) .* sal_triangle,[1 3]),nT,1);
    
    % total energy
    energy = sum(extra.perTriangleRGBError,'all');
end