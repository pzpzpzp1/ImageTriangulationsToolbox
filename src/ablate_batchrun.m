close all; clear all;

boostFactor = 10;
initialHorizontalSamplings = [15 20 25 30 50];
degrees = 0;
maxIters = 500;
saveOut = 1; 
integral1DNsamples = 15;
integral1DNsamplesSubdiv = 50;
edgeSplitResolution = 10;
Nedges2subdivide = 20;
subdivmax = 0; 
subdivisionDamper = 5;
demandedEnergyDensityDrop = 5; windowSize = 20; 
salstrats = SaliencyStrategy.none;
        
foldername = 'images';
files = dir([foldername '/Jcat.jpg']);
for i=1:numel(files)
    [~,name,ext] = fileparts(files(i).name);
    fname = [foldername '/' name ext];
    for j=1:numel(initialHorizontalSamplings)
        initialHorizontalSampling = initialHorizontalSamplings(j);
        degree = 0;
        salstrat = SaliencyStrategy.none;
        

        outputDir = sprintf('output/ablate_%s_init_%d_deg_%d_sal_%d',name,initialHorizontalSampling,degree,0);
        close all;
        ImageTriangulation(fname, ...
            salstrat, boostFactor,...
            initialHorizontalSampling, 0,...
            degree, 0, maxIters, saveOut, outputDir, ... 
            OptStrategy.RMSProp, demandedEnergyDensityDrop, windowSize, ...
            DtStrategy.constrained,...
            integral1DNsamples,...
            integral1DNsamplesSubdiv, edgeSplitResolution, Nedges2subdivide, subdivmax, subdivisionDamper, SubdivisionStrategy.edge, ...
            initialMesh.hexagonal);
    end
end

%%%%%%%%%
boostFactor = 10;
initialHorizontalSamplings = [15];
degrees = 0;
maxIters = 4000;
saveOut = 1; 
integral1DNsamples = 15;
integral1DNsamplesSubdiv = 50;
edgeSplitResolution = 10;
Nedges2subdivides = [5 10 20 30 40];
subdivmax = 10; 
subdivisionDamper = 5;
demandedEnergyDensityDrop = 5; windowSize = 20; 
salstrats = SaliencyStrategy.none;
        
foldername = 'images';
files = dir([foldername '/Jcat.jpg']);
for i=1:numel(files)
    [~,name,ext] = fileparts(files(i).name);
    fname = [foldername '/' name ext];
    for j=1:numel(initialHorizontalSamplings)
        for k=1:numel(Nedges2subdivides)
            Nedges2subdivide = Nedges2subdivides(k);

            initialHorizontalSampling = initialHorizontalSamplings(j);
            degree = 0;
            salstrat = salstrats(1);
            outputDir = sprintf('output/ablate_nsubdivs_%s_init_%d_deg_%d_sal_%d_ndivs_%d',name,initialHorizontalSampling,degree,0,Nedges2subdivide);
            close all;
            ImageTriangulation(fname, ...
                salstrat, boostFactor,...
                initialHorizontalSampling, 0,...
                degree, 0, maxIters, saveOut, outputDir, ... 
                OptStrategy.RMSProp, demandedEnergyDensityDrop, windowSize, ...
                DtStrategy.constrained,...
                integral1DNsamples,...
                integral1DNsamplesSubdiv, edgeSplitResolution, Nedges2subdivide, subdivmax, subdivisionDamper, SubdivisionStrategy.loop, ...
                initialMesh.hexagonal);
        end
    end
end
