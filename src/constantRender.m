function ptc = constantRender(mesh, colors)
    X = mesh.X;
    T = mesh.T;
    ptc = patch('vertices',X ,'faces',T,'edgecolor','none','linewidth',.1,'FaceColor','flat','FaceVertexCdata',uint8(colors),'facealpha',1);
end