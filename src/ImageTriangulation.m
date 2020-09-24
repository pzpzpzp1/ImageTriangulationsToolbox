function [X, T, colors] = ImageTriangulation(fname, ...
        salstrat, boostFactor,...
        initialHorizontalSampling, perturbInit,...
        degree, forceGray, maxIters, saveOut, outputDir, ... 
        optstrat, demandedEnergyDensityDrop, windowSize, ...
        dtstrat,...
        integral1DNsamples,...
        integral1DNsamplesSubdiv, edgeSplitResolution, Nedges2subdivide, subdivmax, subdivisionDamper, substrat, ...
        iMesh...
        )
    
    if nargin == 0
        %% input args
        % INPUT FILE
%         fname = 'images/gradient1.png';
%         fname = 'images/zebra.jpg';
%         fname = 'images/person.jpg';
        fname = 'images/apple.jpg';
%         fname = 'images/face2.jpg';
%         fname = 'images/face3.jpg';
%         fname = 'images/eye.jpg';

        % SALIENCY MAP PARAMETERS
        % salstrat = SaliencyStrategy.none;
        salstrat = SaliencyStrategy.manual;
        % salstrat = SaliencyStrategy.edge; edgediffusion = 5;
        boostFactor = 10;

        % INITIAL MESH DETAILS
        initialHorizontalSampling = 10;
        perturbInit = 0;
        iMesh = initialMesh.hexagonal;
%         iMesh = initialMesh.trim;

        % MISC PARAMETERS
        degree = 0;
        forceGray = 0;
        maxIters = 500;
        saveOut = 0; outputDir = 'output';

        % OPTIMIZATION PARAMETERS
        % optstrat = OptStrategy.none;
        % optstrat = OptStrategy.adaDelta; 
        optstrat = OptStrategy.RMSProp; 
%         demandedEnergyDensityDrop = 0; windowSize = inf;
        demandedEnergyDensityDrop = 5; windowSize = 10; 
        % demandedEnergyDensityDrop = inf; windowSize = 1;
        
        % dtstrat = DtStrategy.none; 
        dtstrat = DtStrategy.constrained;

        % INTEGRATION ACCURACY
        integral1DNsamples = 15;

        % SUBDIVISION PARAMETERS
        integral1DNsamplesSubdiv = 50;
        edgeSplitResolution = 10;
        Nedges2subdivide = 20;
        subdivmax = 10; % times to do subdivision
        subdivisionDamper = 5;
        substrat = SubdivisionStrategy.loop;
        
    end
%% start processing
[path, name, ext] = fileparts(fname);

% load image
img = imread(fname); if size(img,3)==1; img = repmat(img,1,1,3); end
width = size(img,2);
height = size(img,1);
totalArea = width*height;
if forceGray
    greyImg = rgb2gray(img);
    img(:,:,1) = greyImg;
    img(:,:,2) = greyImg;
    img(:,:,3) = greyImg;
end

% get saliency info
salmap = ones(size(img,[1 2]));
if salstrat == SaliencyStrategy.edge
    [Gx1, Gy1] = imgradientxy(img(:,:,1),'prewitt');
    [Gx2, Gy2] = imgradientxy(img(:,:,2),'prewitt');
    [Gx3, Gy3] = imgradientxy(img(:,:,3),'prewitt');
    totalGrad = sqrt(Gx1.^2 + Gy1.^2 + Gx2.^2 + Gy2.^2 + Gx3.^2 + Gy3.^2);
    edgeboost = imgaussfilt(totalGrad, edgediffusion);
    edgeboost = edgeboost/max(max(edgeboost));
    salmap = salmap + edgeboost*boostFactor;
elseif salstrat == SaliencyStrategy.map
    [path,name,ext] = fileparts(fname);
    mapfile = ['saliencymaps/' name ext];
    if exist(mapfile,'file')
        map = double(imread(mapfile));
        map = map/max(max(map));
        salmap = salmap + map*boostFactor;
    else
        display(['No saliency map found at ' mapfile '. Using no saliency to continue.']);
    end
elseif salstrat == SaliencyStrategy.manual
    [path,name,ext] = fileparts(fname);
    mapfile = ['manualsaliency/' name ext];
    if exist(mapfile,'file')
        map = double(imread(mapfile));
        map = map/max(max(map));
        salmap = salmap + map*boostFactor;
    else
        boost = getBoostMap(img);
        imwrite(boost,mapfile);
        if norm(boost,'fro')~=0; boost = boost ./ max(boost(:)); end;
        salmap = salmap + boost*boostFactor;
    end
elseif salstrat == SaliencyStrategy.none
end

% initialize triangulation
if iMesh == initialMesh.hexagonal
    [X,T] = initialHexLatticeMesh(width, height, initialHorizontalSampling);
elseif iMesh == initialMesh.grid
    [X,T] = initialGridMesh(width, height, initialHorizontalSampling, 1);
elseif iMesh == initialMesh.trim
    mappedDensity = (initialHorizontalSampling*.004-.01).^2*9; % maps sampling to density in range [5 - 25].
    [X, T, ~] = imtriangulate(img, mappedDensity); 
end
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
approx0 = Approximator(0);
[~, energy, colors] = approx.computeEnergy(img, mesh, integral1DNsamples,salmap);
areaEnergy0 = getAreaEnergy(mesh, salmap);
areafactor = abs(energy/(areaEnergy0*2));

render(img,mesh,colors,approx,[],salmap);
if saveOut
    outpath = [outputDir '/']; 
    if exist(outpath,'dir')
        display(['SKIPPED: ' outpath]);
        return;
    end
    mkdir(outpath);
    v = VideoWriter([outpath 'mov.avi'],'Motion JPEG AVI'); v.Quality = 80; open(v);
end

%% simulation loop
dt = .4; subdivcount = 1;
subdiviters = zeros(subdivmax,1);
dts = zeros(maxIters,1);
energy = zeros(maxIters,1);
gradnorms = zeros(maxIters,1);
for i=1:maxIters
    %% save output
    title(sprintf('(iter:%d) (subdiviters:%d)', i, subdivcount-1)); drawnow;
    if saveOut
        Xs{i} = X; Ts{i} = T;
        A = getframe(gcf);
        writeVideo(v, A.cdata);
    end

    % update mesh
    mesh = MeshFromXT(X,T);

    %% compute new colors for updated mesh and display
    [extra, approxEnergy, colors, grad] = approx.computeEnergy(img, mesh, integral1DNsamples, salmap);
    [areaEnergy, areaGradient] = getAreaEnergy(mesh, salmap);
    energy(i) = areaEnergy*areafactor + approxEnergy;
    totalGrad = areaGradient*areafactor + grad;
%     totalGrad = areaGradient;
    gradnorms(i) = norm(totalGrad,'fro');

    %% obtain descent direction
    if optstrat==OptStrategy.nonlinearCG
        [descDir, beta] = getNonlinCGDescDir(totalGrad, nonlinearCG);
    elseif optstrat==OptStrategy.adaDelta
        [descDir, rates] = getAdadeltaDescDir(totalGrad, 1); 
    elseif optstrat==OptStrategy.RMSProp
        [descDir, rates] = getAdadeltaDescDir(totalGrad, 2); 
    elseif optstrat==OptStrategy.none
        descDir = -totalGrad;
    end

    render(img,mesh,extra.colorsAlt,approx,[],salmap);

    %% check convergence and either subdivide or stop
    % if energy hasn't dropped significantly since window iterations ago, then energy is 'flat'
    if i > windowSize && energy(i-windowSize) - energy(i) < demandedEnergyDensityDrop*totalArea 
        if subdivcount <= subdivmax
            %% handle bad sliver triangles
            badTriInds = getTrisToCollapse(mesh,extra.perTriangleRGBError,img);
            if numel(badTriInds)~=0
                %{
                figure; image(img); hold all; axis equal; 
                patch('faces',T,'vertices',X,'facecolor','green','facealpha',.1)
                patch('faces',T(badTriInds,:),'vertices',X,'facecolor','none','facealpha',1,'edgecolor','red','linewidth',2)
                %}
                [mesh, reducedInds] = collapseSliverTriangles(mesh, badTriInds);
                if numel(reducedInds)~=0
                    display('sliver collapse created inversion. retrying.');
                    [mesh, furtherreducedInds] = collapseSliverTriangles(mesh, reducedInds);
                    if numel(furtherreducedInds)~=0
                        warning('sliver collapse created inversion. retry failed. skipping collapse.');
                    end
                end
            end        

            %% do subdivision
            subdiviters(subdivcount)=i;
            subdivcount = subdivcount + 1;
            subdivapprox = approx0;
            if substrat == SubdivisionStrategy.edge
                % split mesh via edge cutting
                score = getEdgeSplitScore(mesh, img, subdivapprox, integral1DNsamplesSubdiv, salmap);
                edgeInds = drawEdgesToSplit(Nedges2subdivide, score);
                [X,T] = subdivideMeshEdges(mesh, edgeInds, img, edgeSplitResolution);
            elseif substrat == SubdivisionStrategy.loop
                % split triangles via localized loop subdiv
                score = getTriSplitScore(mesh, img, subdivapprox, integral1DNsamplesSubdiv, salmap);
                [~,perm]=sort(score,'desc'); triInds = perm(1:min(Nedges2subdivide,mesh.nT));
                [X, T] = subdivideMeshTriangles(mesh, triInds);                    
            end
            demandedEnergyDensityDrop = demandedEnergyDensityDrop / subdivisionDamper;
            continue;
        else
            display('Optimization finished!');
            break;
        end
    end

    %% find a good dt to use
    if dtstrat == DtStrategy.constrained
        dt = 1/max(vecnorm(descDir,2,2));
        while any(getTriangleAreas(X+dt*descDir,T)<0)
            dt = dt / 2;
        end
    elseif dtstrat == DtStrategy.onepix
        dt = 1/max(vecnorm(descDir,2,2));
    elseif dtstrat == DtStrategy.none
        if i==1; warning('constant dt is VERY not recommended.'); end
    end
    dts(i) = dt;

    %% update state
    X = X + dt * descDir;
    X = clipVerts(X,width,height);
end
    
render(img,mesh,colors,approx,[],[]);
title(sprintf('(iter:%d) (subdiviters:%d)', i, subdivcount-1));
    
if saveOut
    Xs{i} = X; Ts{i} = T;
    A = getframe(gcf);
    writeVideo(v, A.cdata);
    close(v);
    save([outpath 'XTs.mat'],'Xs','Ts')
end

figure; 
subplot(3,1,1); hold all; title('energy'); plot(energy(1:i-1)); 
for j=1:subdivmax; xline(subdiviters(j)); end
subplot(3,1,2); hold all; title('gradnorm'); plot(gradnorms(1:i-1)); 
subplot(3,1,3); hold all; title('dt'); plot(dts(1:i-1)); 



end
