function ptc = polyRender(mesh, extra, polyparams, renderparams)
    % xys = reshape(extra.samplePoints,[],2); cs = reshape(extra.evalColors,[],3);
    % scatter(xys(:,1), xys(:,2), 50, double(uint8(cs))/255, 'filled')
    
    [ws, ~, bcTris] = getBarycentricSamplingWeights(renderparams.renderResolution); n = size(ws,1); 
    if renderparams.renderResolution ~= renderparams.integral1DNsamples
        samplePoints = mesh.triMaps(ws); 
    else
        samplePoints = extra.samplePoints;
    end
    
    % errors if colstrat isn't amenable to renderparams
    if renderparams.renderWithInteriorColors
        extra.colMaps = poly_triangle_map(extra.colorCPs.int);
    else
        extra.colMaps = poly_triangle_map(extra.colorCPs.full);
    end
    
    evalColors = extra.colMaps(ws);
    
    % render by stacking triangles and using one patch call
    nT = mesh.nT;
    nC = size(samplePoints,1);
    xys = reshape(samplePoints,[],2);
    nt = size(bcTris,1);
    bcTrisFlat = reshape(repmat(reshape(bcTris,nt,1,3),1,nT,1) + (((1:nT)-1)*nC),[],3);
    ptc = patch('vertices',xys,...
        'faces',bcTrisFlat,'edgecolor','none','linewidth',.1,...
        'FaceColor','interp','FaceVertexCdata',...
        uint8(reshape(evalColors,[],3)),'facealpha',1);
    
    title(sprintf('interior colors: %d render res: %s',...
        renderparams.renderWithInteriorColors,...
        num2str(renderparams.renderResolution)));
end