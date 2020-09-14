function ptc = renderMeshEdges(mesh)
    ptc=patch('vertices',mesh.X,'faces',mesh.T,'facecolor','none','edgecolor','c','linewidth',.1);
end