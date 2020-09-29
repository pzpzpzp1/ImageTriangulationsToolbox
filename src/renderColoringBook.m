function figs = renderColoringBook(img, mesh, colors)
    Ks = [5 15 50];
    X = mesh.X; T = mesh.T; nX = size(X,1); nT = size(T,1);
    triCenters = (X(T(:,1),:)+X(T(:,2),:)+X(T(:,3),:))/3;
    
    identity = randi(1000);
    
    for i=1:numel(Ks);
        K = Ks(i);
        pfh = figure; pfh.Units = 'normalized' ;pfh.Position = [0 0 1 1];
        clf; set(pfh,'color','w');
        set(gca, 'YDir','reverse');
        hold all; axis equal; axis off;
        renderMeshEdges(mesh,[.5 .5],'k');
        set(gca,'XTickLabel',{},'YTickLAbel',{},'Box','on')
        [idx, C] = kmeans(colors, K);
        newcolors = C(idx,:);
        patch('vertices',X ,'faces',T,'edgecolor','none','linewidth',.1,'FaceColor','flat','FaceVertexCdata',uint8(newcolors),'facealpha',1);
        title(sprintf('[nX:%d]    [nT:%d]    [K:%d]', nX, nT, K));
%         exportgraphics(pfh,['finished' num2str(K) num2str(identity) '.pdf'])
        
        pfh = figure; pfh.Units = 'normalized'; pfh.Position = [0 0 1 1];
        clf; set(pfh,'color','w');
        X = mesh.X; T = mesh.T; nX = size(X,1); nT = size(T,1);
        set(gca, 'YDir','reverse');
        hold all; axis equal; axis off;
        renderMeshEdges(mesh,[.5 .5],'k');
        set(gca,'XTickLabel',{},'YTickLAbel',{},'Box','on')
        title(sprintf('[nX:%d]    [nT:%d]    [K:%d]', nX, nT, K));
        for j=1:K
            relinds = idx==j;
            text(triCenters(relinds,1),triCenters(relinds,2),num2str(j), 'FontSize', 5);
        end
%         exportgraphics(pfh,['coloringbook' num2str(K) num2str(identity) '.pdf'])
        
        pfh = figure; image(uint8(reshape(C,[],1,3)));
%         exportgraphics(pfh,['colors' num2str(K) num2str(identity) '.pdf'])
        
        
    end
end