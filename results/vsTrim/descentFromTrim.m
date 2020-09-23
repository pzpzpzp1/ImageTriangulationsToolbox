close all; clear all;

initialHorizontalSamplings = 10:5:25;
degree = 0;
maxIters = 300;
saveOut = 1; 
integral1DNsamples = 15;
demandedEnergyDensityDrop = 0; windowSize = inf; 
salstrat = SaliencyStrategy.none;
imesh = initialMesh.trim;

% throwaway variables
boostFactor = 5;
integral1DNsamplesSubdiv = 5;
edgeSplitResolution = 5;
Nedges2subdivide = 5;
subdivmax = 5; 
subdivisionDamper = 5;

foldername = 'images';
files = dir([foldername '/*.jpg']);
for i=1:numel(files)
    [~,name,ext] = fileparts(files(i).name);
    fname = [foldername '/' name ext];
    for j=1:numel(initialHorizontalSamplings)
        initialHorizontalSampling = initialHorizontalSamplings(j);
        
        outputDir = sprintf('results/vsTrim/%s_init_%d_descfromtrim',name,initialHorizontalSampling);
        close all;
        ImageTriangulation(fname, ...
            salstrat, boostFactor,...
            initialHorizontalSampling, 0,...
            degree, 0, maxIters, saveOut, outputDir, ... 
            OptStrategy.RMSProp, demandedEnergyDensityDrop, windowSize, ...
            DtStrategy.constrained,...
            integral1DNsamples,...
            integral1DNsamplesSubdiv, edgeSplitResolution, Nedges2subdivide, subdivmax, subdivisionDamper, SubdivisionStrategy.edge, ...
            imesh);
    end
end

