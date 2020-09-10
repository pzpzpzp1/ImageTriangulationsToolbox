function f1 = constantRender(X, T, colors)
    ptc = patch('vertices',X ,'faces',T,'edgecolor','none','linewidth',.1,'FaceColor','flat','FaceVertexCdata',uint8(colors),'facealpha',1);
end