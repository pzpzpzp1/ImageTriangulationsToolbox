clear all; close all;
files = dir('../../images/*.jpg');

for i=1:numel(files)
    srcname = [files(i).folder '/' files(i).name];
    img = imread(srcname); sal = ones(size(img,[1 2]));
    [~,name,ext]=fileparts(files(i).name);
    
    degree = 0; sal = 1;
    innerfolder = ['../../output_salientarea_slivercollapse_loopsplit/' name '_init_15_deg_0_sal_1'];
    
    try
        load([innerfolder '/XTs.mat']);
        nXs = cellfun(@numel,Xs)/2;
        nTs = cellfun(@numel,Ts)/3;
        iinds = find(nXs(2:end)-nXs(1:end-1));
        iinds = unique([iinds numel(nXs)]);
        approx = Approximator(degree);
        for ii=1:numel(iinds)
            j = iinds(ii);
            X = Xs{j}; T = Ts{j};

            mesh = MeshFromXT(X,T);
            [extra] = approx.computeEnergy(img, mesh, 15, sal);
            colors = extra.colorsAlt; 

            f1 = figure; 
            f1.Visible = 'off'; 
            ax1 = gca; axis off;
            hold all; axis equal; set(gca, 'YDir','reverse');
            approx.render(X, T, extra.colorsAlt); f1.Visible = 'on'; 
            xlim([0 size(img,2)]); ylim([0 size(img,1)]);
            f1.Units = 'normalized'; f1.Position = [0 0 1 1];

            outpicname = ['../../results/compareCompression/' name '_' num2str(j) '.jpg'];
            exportgraphics(f1,outpicname,'Resolution',180);

            outpicname = ['../../results/compareCompression/' name '_' num2str(j) '.txt'];
            SaveConstantTriangulationBinary(outpicname, X, T, colors);

            close(f1);
        end
    catch; end
end

