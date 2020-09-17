function f1 = render(img, mesh, colors, approx, grad)
    %% set up figure such that closing it throws an error, 
    % but also you don't have to initialize a figure outside of this function
    persistent pfh;
    if numel(pfh)==0
        f1 = gcf; 
        f1.Name = 'image triangulations interface';
        pfh = f1;
    end
    try
        figure(pfh);
    catch ex
        % this does mean that running render once can cause an error, and the next call won't.
        clear pfh;
        throw ex;
    end
    
    %% draw on figure pfh
    clf; set(pfh,'color','w');
    pfh.Visible = 'off';
    X = mesh.X; T = mesh.T;
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
    
    set(pfh,'Visible','on');
    drawnow;
end