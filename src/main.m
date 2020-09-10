clear all; close all;

% input arg handling
fname = 'images/sunset.png';
fname = 'images/toucan.png';
initialHorizontalSampling = 5;
degree = 0;
integral1DNsamples = 10;

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
render(img,mesh,colors,approx);













