clear all; close all;
%% input args
% INPUT FILE
% fname = 'images/gradient1.png';
% fname = 'images/zebra.jpg';
% fname = 'images/person.jpg';
% fname = 'images/apple.jpg';
% fname = 'images/face2.jpg';
fname = 'images/face3.jpg';

% SALIENCY MAP PARAMETERS
% salstrat = SaliencyStrategy.none;
% salstrat = SaliencyStrategy.map;
salstrat = SaliencyStrategy.manual;
% salstrat = SaliencyStrategy.edge; edgediffusion = 5;
boostFactor = 10;

% INITIAL MESH DETAILS
initialHorizontalSampling = 15;
perturbInit = 0;

% MISC PARAMETERS
degree = 0;
forceGray = 0;
maxIters = 250;

% OPTIMIZATION PARAMETERS
% optstrat = OptStrategy.none;
% optstrat = OptStrategy.nonlinearCG; nonlinearCG = 2;
optstrat = OptStrategy.adaDelta; rateflag = 2;
% demandedEnergyDensityDrop = 5; windowSize = 5;
demandedEnergyDensityDrop = inf; windowSize = 1;

% dtstrat = DtStrategy.none; 
dtstrat = DtStrategy.onepix;
% dtstrat = DtStrategy.linesearch; lineSearchFactor = 4/3; 

% INTEGRATION ACCURACY
integral1DNsamples = 15;

% SUBDIVISION PARAMETERS
integral1DNsamplesSubdiv = 50;
edgeSplitResolution = 5;
Nedges2subdivide = 5;
subdivmax = 20; % times to do subdivision

%% start processing
% load image
img = imread(fname); if size(img,3)==1; img = repmat(img,1,1,3); end
width = size(img,2);
height = size(img,1);
totalArea = width*height;
if forceGray
    greyImg = rgb2gray(img);
    img(:,:,1) = greyImg;
    img(:,:,2) = greyImg;
    img(:,:,3) = greyImg;
end

% get saliency info
salmap = ones(size(img,[1 2]));
if salstrat == SaliencyStrategy.edge
    [Gx1, Gy1] = imgradientxy(img(:,:,1),'prewitt');
    [Gx2, Gy2] = imgradientxy(img(:,:,2),'prewitt');
    [Gx3, Gy3] = imgradientxy(img(:,:,3),'prewitt');
    totalGrad = sqrt(Gx1.^2 + Gy1.^2 + Gx2.^2 + Gy2.^2 + Gx3.^2 + Gy3.^2);
    edgeboost = imgaussfilt(totalGrad, edgediffusion);
    edgeboost = edgeboost/max(max(edgeboost));
    salmap = salmap + edgeboost*boostFactor;
elseif salstrat == SaliencyStrategy.map
    [path,name,ext] = fileparts(fname);
    mapfile = ['saliencymaps/' name ext];
    if exist(mapfile,'file')
        map = double(imread(mapfile));
        map = map/max(max(map));
        salmap = salmap + map*boostFactor;
    else
        display(['No saliency map found at ' mapfile '. Using no saliency to continue.']);
    end
elseif salstrat == SaliencyStrategy.manual
    boost = getBoostMap(img);
    if norm(boost,'fro')~=0; boost = boost ./ max(boost(:)); end;
    salmap = salmap + boost*boostFactor;
elseif salstrat == SaliencyStrategy.none
    
end

% initialize triangulation
% [X,T] = initialGridMesh(width, height, initialHorizontalSampling, 1);
[X,T] = initialHexLatticeMesh(width, height, initialHorizontalSampling);
mesh = MeshFromXT(X,T);

% perturb interior vertices for more randomness
if perturbInit
    X(mesh.isInterior,:) = X(mesh.isInterior,:) + randn(size(X(mesh.isInterior,:)))*mean(sqrt(2*mesh.triAreas))/15;
    mesh = MeshFromXT(X,T);
end

% initialize approximator
% get initial polyart
% display initial state
approx = Approximator(degree);
[extra, energy, colors] = approx.computeEnergy(img, mesh, integral1DNsamples,salmap);
render(img,mesh,colors,approx,[],salmap);

%% simulation loop
try
    dt = .4; subdivcount = 1;
    subdiviters = zeros(subdivmax,1);
    dts = zeros(maxIters,1);
    energy = zeros(maxIters,1);
    gradnorms = zeros(maxIters,1);
    for i=1:maxIters
        mesh = MeshFromXT(X,T);
        
        %% compute new colors for updated mesh and display
        [extra, energy(i), colors, grad] = approx.computeEnergy(img, mesh, integral1DNsamples, salmap);
        gradnorms(i) = norm(grad,'fro');
        
        %% obtain descent direction
        if optstrat==OptStrategy.nonlinearCG
            [descDir, beta] = getNonlinCGDescDir(grad, nonlinearCG);
        elseif optstrat==OptStrategy.adaDelta
            [descDir, rates] = getAdadeltaDescDir(grad, rateflag);
            shortScatter(X, rates, -1);
        elseif optstrat==OptStrategy.none
            descDir = -grad;
        end
        
        render(img,mesh,extra.colorsAlt,approx,descDir,salmap);
        
        %% check convergence and either subdivide or stop
        % if energy hasn't dropped significantly since window iterations ago, then energy is 'flat'
        if i > windowSize && energy(i-windowSize) - energy(i) < demandedEnergyDensityDrop*totalArea 
            if subdivcount <= subdivmax
                % do subdivision
                subdiviters(subdivcount)=i;
                subdivcount = subdivcount + 1;
                score = getEdgeSplitScore(mesh, img, approx, integral1DNsamplesSubdiv);
                
                sal_edges = double(sampleImage(salmap, getEdgeSamplePoints(mesh,integral1DNsamples)));
                edgeSalBoost = sum(sal_edges,2)/integral1DNsamples;
                
                edgeInds = drawEdgesToSplit(Nedges2subdivide, score.*edgeSalBoost);
                [X,T,dividedEdgeInds] = subdivideMeshEdges(mesh, edgeInds, img, edgeSplitResolution);
                demandedEnergyDensityDrop = demandedEnergyDensityDrop / 10;
                continue;
            else
                display('Optimization finished!');
                break;
            end
        end
        
        %% find a good dt to use
        if dtstrat == DtStrategy.linesearch
            % line search. doesn't work well because energy is not differentiable and all integrals are approximated.
            dt = dt * lineSearchFactor^2;
            nextEnergy = inf;
            while nextEnergy > energy(i)
                dt = dt/lineSearchFactor;
                Xls = X + dt * descDir;
                if any(getTriangleAreas(Xls,T)<0); continue; end
                [~, nextEnergy] = approx.computeEnergy(img, MeshFromXT(Xls,T), integral1DNsamples);
            end
        elseif dtstrat == DtStrategy.onepix
            dt = 1/max(vecnorm(descDir,2,2));
        elseif dtstrat == DtStrategy.none
            if i==1; warning('constant dt is VERY not recommended.'); end
        end
        dts(i) = dt;
        
        %% update state
        X = X + dt * descDir;
    end
catch ex
    erStack = ex.stack;
end
render(img,mesh,colors,approx,[],[]);
figure; 
subplot(3,1,1); hold all; title('energy'); plot(energy(1:i-1)); 
for i=1:subdivmax; xline(subdiviters(i)); end
subplot(3,1,2); hold all; title('gradnorm'); plot(gradnorms(1:i-1)); 
subplot(3,1,3); hold all; title('dt'); plot(dts(1:i-1)); 








