function f1 = render(img, mesh, colors, approx, grad)
    X = mesh.X; T = mesh.T;
    f1 = gcf; f1.Name = 'image triangulations interface';
    clf; set(gcf,'color','w');
    
    f1.Visible = 'off';
    % figure('name','image triangulations interface','Visible','off');
    subplot_er(1,2,1); set(gca, 'YDir','reverse');
    image(img); hold all; axis equal; axis off;
    renderMeshEdges(mesh,[.5 .5]);
    set(gca,'XTickLabel',{},'YTickLAbel',{},'Box','on')
    
    subplot_er(1,2,2);
    hold all; axis equal; axis off;
    set(gca, 'YDir','reverse');
    approx.render(X, T, colors);
    scatter(X(mesh.isXvert,1),X(mesh.isXvert,2),'r.');
    scatter(X(mesh.isYvert,1),X(mesh.isYvert,2),'b.');
    xlim([min(X(:,1)) max(X(:,1))]);
    ylim([min(X(:,2)) max(X(:,2))]);
    if numel(grad)~=0; quiver(X(:,1), X(:,2), grad(:,1), grad(:,2),'c'); end;
    set(gca,'XTickLabel',{},'YTickLAbel',{},'Box','on');
    
    set(gcf,'Visible','on');
    drawnow;
end