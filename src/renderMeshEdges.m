function ptc = renderMeshEdges(mesh,shift,edgecol)
    if nargin < 3
        edgecol = 'c';
    end
    X2 = mesh.X + shift;
    ptc=patch('vertices',X2,'faces',mesh.T,'facecolor','none','edgecolor',edgecol,'linewidth',.1);
end