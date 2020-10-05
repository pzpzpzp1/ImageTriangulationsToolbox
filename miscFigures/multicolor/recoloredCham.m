clear all; close all;

load D:\Documents\MATLAB\ImageTriangulations\output_salientarea_slivercollapse_loopsplit\chameleon_init_25_deg_0_sal_1\XTs.mat;
% load D:\Documents\MATLAB\ImageTriangulations\output_salientarea_slivercollapse_loopsplit\chameleon_init_25_deg_0_sal_2\XTs.mat;
img = imread('images/chameleon.jpg');
i=numel(Ts);

mesh = MeshFromXT(Xs{i},Ts{i}); X=mesh.X; T=mesh.T;
approx = Approximator(0);
[extra, energy, colors] = approx.computeEnergy(img, mesh, 15, []);
load('colordistro.mat');
colordistro0 = colordistro;

artimg = imread('cham.PNG');
figure; image(artimg); hold all; axis equal; title('before');

% artimg = imread('cham.PNG');
% [idx, C] = kmeans(reshape(double(artimg),[],3), 500);
% newcolors = C(idx,:);
% newimg = uint8(reshape(newcolors,size(artimg)));
% colordistro = uint8(C(Cidx,:));
% save('colordistro.mat','colordistro');
% figure; image(artimg); hold all; axis equal; title('before');
% figure; image(newimg); hold all; axis equal; title('after');
% [~,Cidx] = sort(vecnorm(C,2,2)); 
% figure; image(reshape(uint8(C(Cidx,:)),[],1,3))

[idx, C] = kmeans(colors, 10);
newcolors = C(idx,:);
% pfh = figure; %pfh.Units = 'normalized' ;pfh.Position = [0 0 1 1];
% clf; set(pfh,'color','w'); set(gca, 'YDir','reverse'); hold all; axis equal; axis off; set(gca,'XTickLabel',{},'YTickLAbel',{},'Box','on')
% patch('vertices',X ,'faces',T,'edgecolor','none','linewidth',.1,'FaceColor','flat','FaceVertexCdata',uint8(newcolors),'facealpha',1);

specialT = false(mesh.nT,1);
load('removedTris.mat')
specialT(removedTris)=1;
keeptris = ~specialT;
% boundaryColors = unique(newcolors(specialT,:),'rows')
% boundaryTris = ismember(newcolors,boundaryColors,'rows');
% patch('vertices',X ,'faces',T(specialT,:),'edgecolor','red','linewidth',.1,'FaceColor','none');
% edgestokeep = vecnorm(newcolors(mesh.edges2triangles(:,1),:)-newcolors(mesh.edges2triangles(:,2),:),2,2)==0;
% g = graph();g=addedge(g,mesh.edges2triangles(edgestokeep,1),mesh.edges2triangles(edgestokeep,2))
% bins = conncomp(g);
% keeptris = ~ismember(bins,unique(bins(specialT)));

% % colordistro = colordistro0(50:475,:);
% % sampleinds = randsample(size(colordistro,1), size(colors,1),'true');
% % sampcols = colordistro(sampleinds,:);
% % shiftcols = colors - mean(colors) + mean(sampcols);
% % % shiftcols = colors;
% % bias = .4;
% % newcolors = uint8(double(sampcols)*bias + (1-bias)*shiftcols);
% % 
% % satfac = 1.8;
% % HSV = rgb2hsv(double(newcolors)/255);
% % HSV(:, 2) = HSV(:, 2) * satfac;
% % HSV(HSV > 1) = 1; 
% % newcolors = uint8(hsv2rgb(HSV)*255);
% % 
% % pfh = figure; %pfh.Units = 'normalized' ;pfh.Position = [0 0 1 1];
% % clf; set(pfh,'color','w'); set(gca, 'YDir','reverse'); hold all; axis equal; axis off; set(gca,'XTickLabel',{},'YTickLAbel',{},'Box','on')
% % patch('vertices',X ,'faces',T(keeptris,:),'edgecolor','none','linewidth',.1,'FaceColor','flat','FaceVertexCdata',uint8(newcolors(keeptris,:)),'facealpha',1);

colordistro = colordistro0(1:475,:);
pfh = figure; pfh.Units = 'normalized' ;pfh.Position = [0 0 1 1];
clf; set(pfh,'color','w'); set(gca, 'YDir','reverse'); hold all; axis equal; axis off; set(gca,'XTickLabel',{},'YTickLAbel',{},'Box','on')
patch('vertices',X ,'faces',T(keeptris,:),'edgecolor','none','linewidth',.1,'FaceColor','flat','FaceVertexCdata',uint8(colors(keeptris,:)),'facealpha',1);
Ks = [2 32 62 200 400 ];
for i=1:numel(Ks);
    K=Ks(i);
    [idxx] = knnsearch(colordistro, colors, 'k', K);
    [~,idxxx] = max(rand(size(idxx))',[],'linear'); 
    idxxT = idxx';
    idxx = idxxT(idxxx);
    % idxx = idxx(:,2);
    newcolors2 = colordistro(idxx,:);
    pfh = figure; pfh.Units = 'normalized' ;pfh.Position = [0 0 1 1];
    clf; set(pfh,'color','w'); set(gca, 'YDir','reverse'); hold all; axis equal; axis off; set(gca,'XTickLabel',{},'YTickLAbel',{},'Box','on')
    patch('vertices',X ,'faces',T(keeptris,:),'edgecolor','none','linewidth',.1,'FaceColor','flat','FaceVertexCdata',uint8(newcolors2(keeptris,:)),'facealpha',1);
%     patch('vertices',X ,'faces',T(keeptris,:),'edgecolor','none','linewidth',.1,'FaceColor','flat','FaceVertexCdata',uint8(colors(keeptris,:)),'facealpha',1);
end

% removedTris = find(~keeptris)
% save('removedTris.mat','removedTris')

% Xp = [948 322.3];
% Xi=find(vecnorm(mesh.X-Xp,2,2)<.1);
% Tvs = find(any(T==Xi,2))'
% 
% Xp1 = [604.2 246.6];
% Xp2 = [637.6 240.6];
% Xi1=find(vecnorm(mesh.X-Xp1,2,2)<.1);
% Xi2=find(vecnorm(mesh.X-Xp2,2,2)<.1);
% Tvs1 = find(any(T==Xi1,2))'
% Tvs2 = find(any(T==Xi2,2))'
% Tvs = intersect(Tvs1,Tvs2)
