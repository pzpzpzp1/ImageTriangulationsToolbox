
% vn will be (n, nT, 3, 6). 3 for number of edges. 6 for number of velocity values.
function vn_dl  = sampleVdotN_dl(mesh, n)
    T = mesh.T;
    nT = size(T,1);
    
    ws = linspace(0,1,n);
    TEN = mesh.triangleEdgeNormals; 
    
    % load cached constant v sampling indexing tensor. 
    vsampler = loadVSamplingMatrix;
    vsamples = vsampler.*ws';

    % size(vsamples): n  1 3 6 2
    % size(TEN):      1 nT 3 1 2
    
    vn = sum(reshape(vsamples,n,1,3,6,2).*reshape(TEN,1,nT,3,1,2),5);
    % size(vn): n nT 3 6
    
    vn_dl = vn .* reshape(mesh.triangleEdgeLengths,1,nT,3)/n;
end