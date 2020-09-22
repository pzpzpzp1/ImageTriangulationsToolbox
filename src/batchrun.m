close all; clear all;

boostFactor = 10;
initialHorizontalSamplings = [10 15 20 25];
degrees = [0 1];
maxIters = 500;
saveOut = 1; 
integral1DNsamples = 15;
integral1DNsamplesSubdiv = 50;
edgeSplitResolution = 10;
Nedges2subdivide = 20;
subdivmax = 10; 
subdivisionDamper = 5;
demandedEnergyDensityDrop = 5; windowSize = 20; 
areafactors = [1e6 1e5]; 
salstrats = [SaliencyStrategy.none, SaliencyStrategy.manual];
        
foldername = 'images';
files = dir([foldername '/*.jpg']);
for i=1:numel(files)
    [~,name,ext] = fileparts(files(i).name);
    fname = [foldername '/' name ext];
    for j=1:numel(initialHorizontalSamplings)
        initialHorizontalSampling = initialHorizontalSamplings(j);
        for k = 1:numel(degrees)
            degree = degrees(k);
            areafactor = areafactors(k);
            for l = 1:numel(salstrats)
                salstrat = salstrats(l);
                outputDir = sprintf('output/%s_init_%d_deg_%d_sal_%d/',name,initialHorizontalSampling,degree,l);
                ImageTriangulation(fname, ...
                    salstrat, boostFactor,...
                    initialHorizontalSampling, 0,...
                    degree, 0, maxIters, saveOut, outputDir, ... 
                    OptStrategy.RMSProp, demandedEnergyDensityDrop, windowSize, areafactor,...
                    DtStrategy.constrained,...
                    integral1DNsamples,...
                    integral1DNsamplesSubdiv, edgeSplitResolution, Nedges2subdivide, subdivmax, subdivisionDamper ...
                    );
            end
        end
    end
end

