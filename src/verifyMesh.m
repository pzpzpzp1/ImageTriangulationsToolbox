function verifyMesh(mesh)
    X = mesh.X;
    T = mesh.T;
    edgeX = mesh.edgeX;
    faceX = mesh.faceX;
    edges = mesh.edges;
    nT = mesh.nT;
    nX = mesh.nX;
    [~, prePerm] = XedgeXfaceXtoControlPoints(mesh);

    %% display various things about this triangle based on prePerm
    i = randi(nT);
    figure; 
    clf; hold all; axis equal; 
    scatter(prePerm.EX(i,1,:),prePerm.EX(i,2,:),'r','filled')
    for j=1:size(prePerm.EX,3)
        text(prePerm.EX(i,1,j), prePerm.EX(i,2,j), ['<-- e' num2str(j)]);
    end
    scatter(prePerm.TX(i,1,:),prePerm.TX(i,2,:),'k','filled')
    for j=1:3
        text(prePerm.TX(i,1,j), prePerm.TX(i,2,j), ['<-- v' num2str(j)]);
    end
    scatter(prePerm.FX(i,1,:),prePerm.FX(i,2,:),'g','filled')
    for j=1:size(prePerm.FX,3)
        text(prePerm.FX(i,1,j), prePerm.FX(i,2,j), ['<-- f' num2str(j)]);
    end
    title('verts ccw. face cw. edges ccw. per edge, cw.');
    
    
end
