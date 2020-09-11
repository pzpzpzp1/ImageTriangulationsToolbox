clear all; close all;

% input arg handling
fname = 'images/toucan.png';
initialHorizontalSampling = 3;
degree = 0;
dt0 = .000001; % initial dt
integral1DNsamples = 30;
maxIters = 1000;

% load image
img = imread(fname);
width = size(img,2);
height = size(img,1);

% initialize triangulation
[X,T] = initialMesh(img, initialHorizontalSampling);
mesh = MeshFromXT(X,T);

% initialize approximator
approx = Approximator(degree);

% get initial polyart
[energy, colors] = constantComputeEnergy(img, mesh, integral1DNsamples);

% display initial state
render(img,mesh,colors,approx,[]);

%% simulation loop
dt = dt0;
for i=1:maxIters
    mesh = MeshFromXT(X,T);
    [energy, colors] = constantComputeEnergy(img, mesh, integral1DNsamples);
    grad = approx.computeGradient(img, mesh, integral1DNsamples);
    
    render(img,mesh,colors,approx,grad);
    
    X = X - dt * grad;
end









