clear all; close all;

%% input args
% fname = 'images/apple.png';
% fname = 'images/sunset.png';
% fname = 'images/circle.png';
fname = 'images/BW.png';
% fname = 'images/toucan.png';
initialHorizontalSampling = 10;
degree = 0;
dt0 = 5e-9; % initial dt
integral1DNsamples = 20;
integral1DNsamplesSubdiv = 50;
edgeSplitResolution = 10;
maxIters = 1;
Nedges2subdivide = 50;
forceGray = 0;
perturbInit = 1;

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
    dt = dt0;
    energy = zeros(maxIters,1);
    gradnorms = zeros(maxIters,1);
    for i=1:maxIters
        mesh = MeshFromXT(X,T);
        
        [extra, energy(i), colors, grad] = approx.computeEnergy(img, mesh, integral1DNsamples);
        gradnorms(i) = norm(grad);
        
        render(img,mesh,extra.colorsAlt,approx,grad);

        dt = 1/max(vecnorm(grad,2,2)); % always at most 1 pixel distance traveled per vertex.
        if norm(grad)==0; break; end;
        X = X - dt * grad;
    end
catch ex
    erStack = ex.stack;
end
render(img,mesh,colors,approx,[]);
figure; 
subplot(1,2,1); hold all; title('energy'); plot(energy(1:i)); 
subplot(1,2,2); hold all; title('gradnorm'); plot(gradnorms(1:i)); 

%% do subdivision method
score = getEdgeSplitScore(mesh, img, approx, integral1DNsamplesSubdiv);
edgeInds = drawEdgesToSplit(Nedges2subdivide, score);
[X,T,dividedEdgeInds] = subdivideMeshEdges(mesh, edgeInds, img, edgeSplitResolution);
mesh = MeshFromXT(X,T);
[extra, energy(i), colors] = approx.computeEnergy(img, mesh, integral1DNsamples);
render(img,mesh,extra.colorsAlt,approx,[]);








