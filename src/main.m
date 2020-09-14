clear all; close all;

% input arg handling
fname = 'images/sunset.png';
% fname = 'images/gradientVerticalGray.png';
initialHorizontalSampling = 20;
degree = 1;
% dt0 = 1e-6; % initial dt
dt0 = 5e-8; % initial dt
integral1DNsamples = 20;
% integral1DNsamples = 500;
maxIters = 1000;
showgrad = 1;

% load image
img = imread(fname);
width = size(img,2);
height = size(img,1);

% initialize triangulation
% [X,T] = initialGridMesh(width, height, initialHorizontalSampling);
[X,T] = initialHexLatticeMesh(width, height, initialHorizontalSampling);
mesh = MeshFromXT(X,T);

% perturb interior vertices for more randomness
% X(mesh.isInterior,:) = X(mesh.isInterior,:) + randn(size(X(mesh.isInterior,:)))*width/(20*initialHorizontalSampling);
% mesh = MeshFromXT(X,T);

% initialize approximator
approx = Approximator(degree);

% get initial polyart
[energy, colors] = approx.computeEnergy(img, mesh, integral1DNsamples);

% display initial state
render(img,mesh,colors,approx,[]);

%% simulation loop
dt = dt0;
energy = zeros(maxIters,1);
for i=1:maxIters
    mesh = MeshFromXT(X,T);
    [energy(i), colors, grad] = approx.computeEnergy(img, mesh, integral1DNsamples);
    
    if showgrad;
        render(img,mesh,colors,approx,grad);
    else
        render(img,mesh,colors,approx,[]);
    end
    
    dt = 1/max(vecnorm(grad,2,2)); % always at most 1 pixel distance traveled per vertex.
    X = X - dt * grad;
end
render(img,mesh,colors,approx,[]);
figure; plot(energy)







