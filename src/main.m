clear all; close all;

%% input args
% fname = 'images/appleGray.PNG';
% fname = 'images/apple.png';
% fname = 'images/sunset.png';
% fname = 'images/circle.png';
% fname = 'images/BW.png';
% fname = 'images/s2by3.PNG';
fname = 'images/toucan.png';
% fname = 'images/gradientVerticalGray.png';
% fname = 'images/gradientDiagGray.png';
% fname = 'images/gradientHorizontalGray.png';
initialHorizontalSampling = 25;
degree = 1;
dt0 = 5e-9; % initial dt
integral1DNsamples = 15;
maxIters = 1000;
showgrad = 1;
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
approx = Approximator(degree);

% get initial polyart
[energy, colors] = approx.computeEnergy(img, mesh, integral1DNsamples);

% display initial state
render(img,mesh,colors,approx,[]);

%% simulation loop
try
    dt = dt0;
    energy = zeros(maxIters,1);
    gradnorms = zeros(maxIters,1);
    for i=1:maxIters
        mesh = MeshFromXT(X,T);
        [energy(i), colors, grad] = approx.computeEnergy(img, mesh, integral1DNsamples);
        gradnorms(i) = norm(grad);
        
        if showgrad;
            render(img,mesh,colors,approx,grad);
        else
            render(img,mesh,colors,approx,[]);
        end

        dt = 1/max(vecnorm(grad,2,2)); % always at most 1 pixel distance traveled per vertex.
        if norm(grad)==0; break; end;
        X = X - dt * grad;
    end
catch ex
    erStack = ex.stack;
end
render(img,mesh,colors,approx,[]);
figure; plot(energy(1:i)); title('energy');
figure; plot(gradnorms(1:i)); title('gradnorm');







