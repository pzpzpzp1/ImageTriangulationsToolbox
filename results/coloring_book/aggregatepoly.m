clear all; close all;

% load D:\Documents\MATLAB\ImageTriangulations\output_salientarea_slivercollapse_loopsplit\fish_init_15_deg_0_sal_2\XTs.mat;
% img = imread('images/fish.jpg');
% % i=217;
% i=127;
% Ks = [10 100];

% load D:\Documents\MATLAB\ImageTriangulations\output_salientarea_slivercollapse_loopsplit\geyser_init_15_deg_0_sal_1\XTs.mat;
% img = imread('images/geyser.jpg');
% i=151;
% Ks = [10 100];


load D:\Documents\MATLAB\ImageTriangulations\output_salientarea_slivercollapse_loopsplit\cupcake_init_15_deg_0_sal_2\XTs.mat;
img = imread('images/cupcake.jpg');
i=150;
% i=115;
K = 100;

mesh = MeshFromXT(Xs{i},Ts{i}); X=mesh.X; T=mesh.T;
approx = Approximator(0);
[extra, energy, colors] = approx.computeEnergy(img, mesh, 15, []);

% pfh = figure; pfh.Units = 'normalized' ;pfh.Position = [0 0 1 1];
% clf; set(pfh,'color','w');set(gca, 'YDir','reverse');hold all; axis equal; axis off;set(gca,'XTickLabel',{},'YTickLAbel',{},'Box','on')
% renderMeshEdges(mesh,[.5 .5],'k');

[idx, C] = kmeans(colors, K);
newcolors = C(idx,:);

edgestokeep = vecnorm(newcolors(mesh.edges2triangles(:,1),:)-newcolors(mesh.edges2triangles(:,2),:),2,2)==0;
g = graph();g=addedge(g,mesh.edges2triangles(edgestokeep,1),mesh.edges2triangles(edgestokeep,2))
bins = conncomp(g);
uniquebins = unique(bins);

singletris = [];
doubletris = [];
tripletris = [];
aggedpoly = [];
for i=1:numel(uniquebins)
    bind = uniquebins(i);
    tristojoin = find(bins == bind);
    aggedpoly(i) = numel(tristojoin);
    if aggedpoly(i)==1
        singletris(end+1) = tristojoin;
    elseif aggedpoly(i)==2
        doubletris(end+1,:) = tristojoin;
    elseif aggedpoly(i)==3
        tripletris(end+1,:) = tristojoin;
    end
end
unique(aggedpoly)

singletris = setdiff(1:mesh.nT,unique([doubletris(:); tripletris(:)]));

ec = 'c'; lw = 1;
pfh = figure; pfh.Units = 'normalized' ;pfh.Position = [0 0 1 1]; clf; set(pfh,'color','w');set(gca, 'YDir','reverse');hold all; axis equal; axis off;set(gca,'XTickLabel',{},'YTickLAbel',{},'Box','on')
patch('vertices',X ,'faces',T,'edgecolor','none','linewidth',.1,'FaceColor','flat','FaceVertexCdata',uint8(newcolors),'facealpha',1);
patch('vertices',X ,'faces',T(singletris,:),'edgecolor',ec,'linewidth',lw,'FaceColor','none');
for i=1:size(doubletris,1)
    twoTris = doubletris(i,:);
    alledges = reshape(mesh.triangles2edges(twoTris,:),[],1);
    uniqueedges = unique(alledges);
      
    [~,col] = ismember(uniqueedges,alledges);
    repedges = alledges;
    repedges(col) = [];
    singleedges = setdiff(alledges,repedges);
    patch('vertices',X ,'faces',mesh.edges(singleedges,[1 2 1]),'edgecolor',ec,'linewidth',lw,'FaceColor','none');
end
for i=1:size(tripletris,1)
    twoTris = tripletris(i,:);
    alledges = reshape(mesh.triangles2edges(twoTris,:),[],1);
    uniqueedges = unique(alledges);
      
    [~,col] = ismember(uniqueedges,alledges);
    repedges = alledges;
    repedges(col) = [];
    singleedges = setdiff(alledges,repedges);
    patch('vertices',X ,'faces',mesh.edges(singleedges,[1 2 1]),'edgecolor',ec,'linewidth',lw,'FaceColor','none');
end

fa = .4;
pfh = figure; pfh.Units = 'normalized' ;pfh.Position = [0 0 1 1]; clf; set(pfh,'color','w');set(gca, 'YDir','reverse');hold all; axis equal; axis off;set(gca,'XTickLabel',{},'YTickLAbel',{},'Box','on'); 
title('polygons highlighted');
patch('vertices',X ,'faces',T,'edgecolor','none','linewidth',.1,'FaceColor','flat','FaceVertexCdata',uint8(newcolors),'facealpha',1);
patch('vertices',X ,'faces',T(doubletris,:),'edgecolor','none','linewidth',.1,'FaceColor','red','facealpha',fa);
patch('vertices',X ,'faces',T(tripletris,:),'edgecolor','none','linewidth',.1,'FaceColor','blue','facealpha',fa);
for i=1:size(doubletris,1)
    twoTris = doubletris(i,:);
    alledges = reshape(mesh.triangles2edges(twoTris,:),[],1);
    uniqueedges = unique(alledges);
      
    [~,col] = ismember(uniqueedges,alledges);
    repedges = alledges;
    repedges(col) = [];
    singleedges = setdiff(alledges,repedges);
    patch('vertices',X ,'faces',mesh.edges(singleedges,[1 2 1]),'edgecolor','k','linewidth',lw,'FaceColor','none');
end
for i=1:size(tripletris,1)
    twoTris = tripletris(i,:);
    alledges = reshape(mesh.triangles2edges(twoTris,:),[],1);
    uniqueedges = unique(alledges);
      
    [~,col] = ismember(uniqueedges,alledges);
    repedges = alledges;
    repedges(col) = [];
    singleedges = setdiff(alledges,repedges);
    patch('vertices',X ,'faces',mesh.edges(singleedges,[1 2 1]),'edgecolor','k','linewidth',lw,'FaceColor','none');
end


pfh = figure; pfh.Units = 'normalized' ;pfh.Position = [0 0 1 1]; clf; set(pfh,'color','w');set(gca, 'YDir','reverse');hold all; axis equal; axis off;set(gca,'XTickLabel',{},'YTickLAbel',{},'Box','on'); 
title('full');
patch('vertices',X ,'faces',T,'edgecolor','none','linewidth',.1,'FaceColor','flat','FaceVertexCdata',uint8(newcolors),'facealpha',1);


