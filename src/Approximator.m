% This used to divert to different functions for different parameters, but now it's unified in the poly functions.
% basically just a reduction in parameters passed around, but could be refactored away.
function approx = Approximator(polyparams, colstrat)
    % poly triangulation, varied color degree
    approx.computeEnergy = @(img, mesh, integral1DNsamples, salmap) polyComputeEnergy(img, mesh, integral1DNsamples, salmap, polyparams, colstrat);
    approx.render = @(mesh,colors,rp) polyRender(mesh,colors,polyparams,rp);
end