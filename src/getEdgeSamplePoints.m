% edgeSamplePoints: [nT, n, 3(edges), 2(xy)]
function edgeSamplePoints = getEdgeSamplePoints(mesh, n)
    X = mesh.X; T=mesh.T; nT=size(T,1);
    
    edgeWs = linspace(0,1,n);
    edgeSamplePoints = zeros(nT,n,3,2);
    v1 = reshape(X(T(:,1),:),nT,1,1,2);
    v2 = reshape(X(T(:,2),:),nT,1,1,2);
    v3 = reshape(X(T(:,3),:),nT,1,1,2);
    edgeSamplePoints(:,:,1,:) = v1.*(1-edgeWs) + v2.*edgeWs;
    edgeSamplePoints(:,:,2,:) = v2.*(1-edgeWs) + v3.*edgeWs;
    edgeSamplePoints(:,:,3,:) = v3.*(1-edgeWs) + v1.*edgeWs;
end