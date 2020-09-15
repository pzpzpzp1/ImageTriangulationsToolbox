function ptc = renderMeshEdges(mesh,shift)
    X2 = mesh.X + shift;
    ptc=patch('vertices',X2,'faces',mesh.T,'facecolor','none','edgecolor','c','linewidth',.1);
end