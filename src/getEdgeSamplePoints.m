% edgeSamplePoints: [nT, n, 3(edges), 2(xy)]
function edgeSamplePoints = getEdgeSamplePoints(mesh, n)
    X = mesh.X; T=mesh.T; nT=size(T,1);
    
    edgeWs = linspace(0,1,n);
    v123 = reshape(X(T(:),:),nT,1,3,2);
    v231 = reshape(X(T(:,[2 3 1]),:),nT,1,3,2);
    edgeSamplePoints = v123.*(1-edgeWs) + v231.*edgeWs;
    
end








