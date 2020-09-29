function f1 = render(img, mesh, colors, approx, grad, salmap)
    %% set up figure such that closing it throws an error, 
    % but also you don't have to initialize a figure outside of this function
    persistent pfh;
    if numel(pfh)==0 || ~isvalid(pfh)
        f1 = gcf; 
        f1.Units = 'normalized';
        f1.OuterPosition = [0 0 1 1];
        f1.Name = 'image triangulations interface';
        pfh = f1;
    end
    
    %% draw on figure pfh
    clf; set(pfh,'color','w');
    pfh.Visible = 'off';
    X = mesh.X; T = mesh.T;
    ax1 = subplot_er(1,2,1); set(gca, 'YDir','reverse');
    image(img); hold all; axis equal; axis off;
    renderMeshEdges(mesh,[.5 .5]);
    set(gca,'XTickLabel',{},'YTickLAbel',{},'Box','on')
    if numel(salmap)~=0 && norm(salmap-1,'fro')~=0
        boost = (salmap-1); imh = image(boost./max(boost(:))*256); imh.AlphaData = .3*boost/max(boost(:));
    end
    title(sprintf('(X:%d) (T:%d)',mesh.nX, mesh.nT));
    
    
    ax2 = subplot_er(1,2,2);
    hold all; axis equal; axis off;
    set(gca, 'YDir','reverse');
    approx.render(X, T, colors);
%     scatter(X(mesh.isXvert,1),X(mesh.isXvert,2),'r.');
%     scatter(X(mesh.isYvert,1),X(mesh.isYvert,2),'b.');
    xlim([min(X(:,1)) max(X(:,1))]);
    ylim([min(X(:,2)) max(X(:,2))]);
    set(gca,'XTickLabel',{},'YTickLAbel',{},'Box','on')
    if numel(grad)~=0; quiver(X(:,1), X(:,2), grad(:,1), grad(:,2),'c'); end;
        
    set(pfh,'Visible','on');
    hlink = linkprop([ax1,ax2],{'XLim','YLim','ZLim'}); 
    pfh.UserData = hlink;
    drawnow;
end