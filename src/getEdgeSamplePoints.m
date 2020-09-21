% edgeSamplePoints: [nE, n, 2(xy)]
function edgeSamplePoints = getEdgeSamplePoints(mesh, n)
    X = mesh.X; T=mesh.T; nT=size(T,1);
    E = mesh.edges; nE = size(E,1);
    
    edgeWs = linspace(0,1,n);
    v1 = reshape(X(E(:,1),:),nE,1,2);
    v2 = reshape(X(E(:,1),:),nE,1,2);
    
    edgeSamplePoints = v1.*(1-edgeWs) + v2.*edgeWs;
    
end








