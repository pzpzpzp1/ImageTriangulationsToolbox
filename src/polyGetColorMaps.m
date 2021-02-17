function colormaps = polyGetColorMaps(img, mesh, integral1DNsamples, salmap, polyparams)
    if numel(salmap)==0
        salmap = img(:,:,1)*0+1;
    end
    
    X = mesh.X; T = mesh.T; nT = size(T,1); 
    % generate sample locations in barycentric coords
    [ws, interiorInds, bcTris] = getBarycentricSamplingWeights(integral1DNsamples); n = size(ws,1); 
    [samplePoints, ~, geodets] = mesh.triMaps(ws); % world space points, and jacobian det
    
    % perform sampling
    f_triangle = double(sampleImage(img, samplePoints)); % image colors
    sal_triangle = double(sampleImage(salmap, samplePoints)); % saliency map values
    
    % get colors per triangle
    [extra.colorCPsAlt, colorCPs] = polyGetColorsFromSamples(ws, f_triangle, sal_triangle, interiorInds, polyparams, geodets);
    extra.colMaps = poly_triangle_map(colorCPs);
    extra.colorCPs = colorCPs;
    extra.samplePoints = samplePoints;
    extra.evalColors = extra.colMaps(ws);
    extra.bcTris = bcTris;
    
    % compute energy integrated over triangles
    extra.perTriangleRGBError = reshape(sum((extra.evalColors - f_triangle).^2 .* abs(geodets) .* sal_triangle,[1 3]),nT,1);
    energy = sum(extra.perTriangleRGBError,'all');
    
end