clear all; close all;

% load D:\Documents\MATLAB\ImageTriangulations\output_salientarea_slivercollapse_loopsplit\chameleon_init_25_deg_0_sal_1\XTs.mat;
% img = imread('images/chameleon.jpg');

load D:\Documents\MATLAB\ImageTriangulations\output_salientarea_slivercollapse_loopsplit\turkey_init_15_deg_0_sal_2\XTs.mat;
img = imread('images/turkey.jpg');

i=numel(Xs);
mesh = MeshFromXT(Xs{i},Ts{i});
approx = Approximator(0);
[extra, energy, colors] = approx.computeEnergy(img, mesh, 15, []);

renderColoringBook(img, mesh, extra.colorsAlt);





