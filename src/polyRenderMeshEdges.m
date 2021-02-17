function ptc = polyRenderMeshEdges(mesh,shift)
    [ws, ~, bcTris] = getBarycentricSamplingWeights(mesh.gdeg+1); n = size(ws,1); 
    samplePoints = mesh.triMaps(ws);
    xys = reshape(samplePoints,[],2);
    nT = mesh.nT;    nC = size(samplePoints,1);    nt = size(bcTris,1);
    bcTrisFlat = reshape(repmat(reshape(bcTris,nt,1,3),1,nT,1) + (((1:nT)-1)*nC),[],3);
    ptc = patch('vertices',xys + shift,...
        'faces',bcTrisFlat,'edgecolor','c','linewidth',.1,...
        'FaceColor','none');
    
    X2 = mesh.X + shift;
    edgeX2 = reshape(permute(mesh.edgeX + shift,[1 3 2]),[],2);
    faceX2 = reshape(permute(mesh.faceX + shift,[1 3 2]),[],2);
    % ptc = patch('vertices',X2,'faces',mesh.T,'facecolor','none','edgecolor',edgecol,'linewidth',.1);
    scatter(X2(:,1),X2(:,2),20,'k','filled')
    scatter(X2(:,1),X2(:,2),15,'r','filled')
    scatter(edgeX2(:,1),edgeX2(:,2),20,'k','filled')
    scatter(edgeX2(:,1),edgeX2(:,2),15,'b','filled')
    scatter(faceX2(:,1),faceX2(:,2),20,'k','filled')
    scatter(faceX2(:,1),faceX2(:,2),15,'g','filled')
end