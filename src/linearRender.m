function ptc = linearRender(X, T, colors)
    colorslin = uint8(reshape(colors,[],3));
    newX = X(reshape(T,[],1),:);
    newT = reshape(1:size(newX,1),[],3);
    ptc = patch('vertices', newX ,'faces',newT,'edgecolor','none','linewidth',.1,...
        'FaceColor','interp','FaceVertexCdata',colorslin,'facealpha',1);
end