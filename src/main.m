clear all; close all;
%% input args
% INPUT FILE
% fname = 'images/apple.png';
% fname = 'images/sunset.png';
% fname = 'images/circle.png';
% fname = 'images/BW.png';
fname = 'images/toucan.png';

% INITIAL MESH DETAILS
initialHorizontalSampling = 30;
perturbInit = 1;

% MISC PARAMETERS
degree = 0;
forceGray = 0;
maxIters = 500;

% OPTIMIZATION PARAMETERS
% optstrat = OptStrategy.none;
% optstrat = OptStrategy.nonlinearCG; nonlinearCG = 2;
optstrat = OptStrategy.adaDelta; rateflag = 2;

% INTEGRATION ACCURACY
integral1DNsamples = 20;

% SUBDIVISION PARAMETERS
integral1DNsamplesSubdiv = 50;
edgeSplitResolution = 10;
Nedges2subdivide = 50;


%% start processing
% load image
img = imread(fname); 
width = size(img,2);
height = size(img,1);

if forceGray
    greyImg = rgb2gray(img);
    img(:,:,1) = greyImg;
    img(:,:,2) = greyImg;
    img(:,:,3) = greyImg;
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
[extra, energy, colors] = approx.computeEnergy(img, mesh, integral1DNsamples);
render(img,mesh,colors,approx,[]);

%% simulation loop
try
    energy = zeros(maxIters,1);
    gradnorms = zeros(maxIters,1);
    for i=1:maxIters
        mesh = MeshFromXT(X,T);
        
        % compute new colors for updated mesh and display
        [extra, energy(i), colors, grad] = approx.computeEnergy(img, mesh, integral1DNsamples);
        gradnorms(i) = norm(grad,'fro');
        
        % obtain descent direction
        if optstrat==OptStrategy.nonlinearCG
            [descDir, beta] = getNonlinCGDescDir(grad, nonlinearCG);
        elseif optstrat==OptStrategy.adaDelta
            [descDir, rates] = getAdadeltaDescDir(grad, rateflag);
            FancyScatter(X, rates, -1);
        elseif optstrat==OptStrategy.none
            descDir = -grad;
        end
        
        render(img,mesh,extra.colorsAlt,approx,descDir);
        
        % todo: line search
        dt = 1/max(vecnorm(descDir,2,2)); % always at most 1 pixel distance traveled per vertex.
        if norm(grad)==0; break; end;
        X = X + dt * descDir;
    end
catch ex
    erStack = ex.stack;
end
% render(img,mesh,colors,approx,[]);
figure; 
subplot(1,2,1); hold all; title('energy'); plot(energy(1:i-1)); 
subplot(1,2,2); hold all; title('gradnorm'); plot(gradnorms(1:i-1)); 

%% do subdivision method
%{
score = getEdgeSplitScore(mesh, img, approx, integral1DNsamplesSubdiv);
edgeInds = drawEdgesToSplit(Nedges2subdivide, score);
[X,T,dividedEdgeInds] = subdivideMeshEdges(mesh, edgeInds, img, edgeSplitResolution);
mesh = MeshFromXT(X,T);
[extra, energy(i), colors] = approx.computeEnergy(img, mesh, integral1DNsamples);
render(img,mesh,extra.colorsAlt,approx,[]);
%}







