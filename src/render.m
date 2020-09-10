function f1 = render(img, mesh, colors, approx)
    X = mesh.X; T = mesh.T;
    
    figure('name','image triangulations interface','Visible','off');
    subplot_er(1,2,1); set(gca, 'YDir','reverse')
    image(img); hold all; axis equal; axis off;
    set(gca,'XTickLabel',{},'YTickLAbel',{},'Box','on')
    
    subplot_er(1,2,2);
    hold all; axis equal; axis off;
    set(gca, 'YDir','reverse')
    approx.render(X, T, colors)
    set(gca,'XTickLabel',{},'YTickLAbel',{},'Box','on')
    
    set(gcf,'Visible','on')
end