clear all; close all;

% load D:\Documents\MATLAB\ImageTriangulations\output_salientarea_slivercollapse_loopsplit\cupcake_init_15_deg_0_sal_2\XTs.mat;
% img = imread('images/cupcake.jpg');
% i=115;
% Ks = [10 100];

load D:\Documents\MATLAB\ImageTriangulations\output_salientarea_slivercollapse_loopsplit\fish_init_15_deg_0_sal_2\XTs.mat;
img = imread('images/fish.jpg');
% i=217;
i=127;
Ks = [10 100];

load D:\Documents\MATLAB\ImageTriangulations\output_salientarea_slivercollapse_loopsplit\geyser_init_15_deg_0_sal_1\XTs.mat;
img = imread('images/geyser.jpg');
i=151;
Ks = [10 100];


mesh = MeshFromXT(Xs{i},Ts{i});
approx = Approximator(0);
[extra, energy, colors] = approx.computeEnergy(img, mesh, 15, []);

renderColoringBook(img, mesh, extra.colorsAlt, Ks);





