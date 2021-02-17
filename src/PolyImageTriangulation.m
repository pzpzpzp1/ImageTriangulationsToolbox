function [X, T, colors, timedata] = PolyImageTriangulation(fname, ...
        salstrat, boostFactor,...
        initialHorizontalSampling, perturbInit,...
        degree, forceGray, maxIters, saveOut, outputDir, ... 
        optstrat, demandedEnergyDensityDrop, windowSize, ...
        dtstrat,...
        integral1DNsamples,...
        integral1DNsamplesSubdiv, edgeSplitResolution, Nedges2subdivide, subdivmax, subdivisionDamper, substrat, ...
        iMesh...
        )
    timedata.energycomptime = [];
    timedata.nXs = [];
    timedata.nTs = [];
    timedata.subdivtime = [];
    functionstarttime = tic;
    if nargin == 0
        % close all;
        
        %% input args
        % INPUT FILE
%         fname = 'images/gradient1.png';
        fname = 'images/redball.jpg';
%         fname = 'images/zebra.jpg';
%         fname = 'images/person.jpg';
%         fname = 'images/Jcat.jpg';
%         fname = 'images/Jflower.jpg';
%         fname = 'images/fish.jpg';
%         fname = 'images/geyser.jpg';
%         fname = 'images/apple.jpg';
%         fname = 'images/checkerboard.JPG';
%         fname = 'images/face2.jpg';
%         fname = 'images/face3.jpg';
%         fname = 'images/eye.jpg';

        % SALIENCY MAP PARAMETERS
        salstrat = SaliencyStrategy.none;
%         salstrat = SaliencyStrategy.manual;
        % salstrat = SaliencyStrategy.edge; edgediffusion = 5;
        boostFactor = 10;

        % INITIAL MESH DETAILS
        initialHorizontalSampling = 7;
        perturbInit = 1;
        iMesh = initialMesh.hexagonal;
%         iMesh = initialMesh.trim;

        % MISC PARAMETERS
        forceGray = 0;
        maxIters = 500;
        saveOut = 0; outputDir = 'outputcache';
        polyparams.gdeg = 1;
        polyparams.cdeg = 0; % best as 0 1 or 2.
        
        % OPTIMIZATION PARAMETERS
%         optstrat = OptStrategy.none;
        % optstrat = OptStrategy.adaDelta; 
        optstrat = OptStrategy.RMSProp; 
%         demandedEnergyDensityDrop = 0; windowSize = inf;
        demandedEnergyDensityDrop = 5; windowSize = 10; 
%         demandedEnergyDensityDrop = inf; windowSize = 1;
        
        % dtstrat = DtStrategy.none; 
        % dtstrat = DtStrategy.constrained;
        dtstrat = DtStrategy.onepix;

        % INTEGRATION ACCURACY
        integral1DNsamples = 10;
        
        % RENDER PARAMETERS
        renderparams.renderWithInteriorColors = false;
        renderparams.renderResolution = integral1DNsamples;
        renderparams.integral1DNsamples = integral1DNsamples;

        % SUBDIVISION PARAMETERS
        integral1DNsamplesSubdiv = 50;
        edgeSplitResolution = 10;
        Nedges2subdivide = 20;
        subdivmax = 10; % times to do subdivision
        subdivisionDamper = 5;
        substrat = SubdivisionStrategy.loop;
%         substrat = SubdivisionStrategy.edge;
        
        colstrat = colorStrategy.full;
        % colstrat = colorStrategy.int;
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
mesh = PolyMeshFromXT(X,T,[],[],[],false,polyparams.gdeg);
edgeX = mesh.edgeX;
faceX = mesh.faceX;
    
% perturb interior vertices for more randomness
if perturbInit
    randscale = min(width,height)/integral1DNsamples/60;
    X(mesh.isInterior,:) = X(mesh.isInterior,:) + randn(size(X(mesh.isInterior,:)))*randscale;
    faceX = faceX + randn(size(faceX))*randscale;
    edgeX(~mesh.isBoundaryEdge,:,:) = edgeX(~mesh.isBoundaryEdge,:,:) + randn(size(edgeX(~mesh.isBoundaryEdge,:,:)))*randscale;
    mesh = PolyMeshFromXT(X,T,mesh.edges,edgeX,faceX,false,polyparams.gdeg);
end

% initialize approximator
% get initial polyart
% display initial state
approx = Approximator(polyparams, colstrat);
[extra, energy0] = approx.computeEnergy(img, mesh, integral1DNsamples, salmap);
renderparams0.renderWithInteriorColors = false;renderparams0.renderResolution = integral1DNsamples;renderparams0.integral1DNsamples = integral1DNsamples;


figure; render(img, mesh, extra, approx, [], salmap, renderparams);
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
    if i~=1
        mesh = PolyMeshFromXT(X,T,mesh.edges,edgeX,faceX,false,polyparams.gdeg);
    end
    
    energycomptimeStart = tic;
    %% compute new colors for updated mesh and display
    [extra, approxEnergy, grad] = approx.computeEnergy(img, mesh, integral1DNsamples, salmap);
    %[areaEnergy, areaGradient] = getAreaEnergy(mesh, salmap);
    %energy(i) = areaEnergy*areafactor + approxEnergy;
    %totalGrad = areaGradient*areafactor + grad;
    
    energy(i) = approxEnergy;
    totalGrad = grad;
    
    gradnorms(i) = norm([totalGrad.faceGrad(:);totalGrad.edgeGrad(:);totalGrad.Xgrad(:)]);
    timedata.energycomptime(i) = toc(energycomptimeStart);
    timedata.nXs(i) = size(X,1);
    timedata.nTs(i) = size(T,1);
    
    %% obtain descent direction
%     if optstrat==OptStrategy.nonlinearCG
%         [descDir, beta] = getNonlinCGDescDir(totalGrad, nonlinearCG);
%     elseif optstrat==OptStrategy.adaDelta
%         [descDir, rates] = getAdadeltaDescDir(totalGrad, 1); 
%     elseif optstrat==OptStrategy.RMSProp
%         [descDir, rates] = getAdadeltaDescDir(totalGrad, 2); 
%     elseif optstrat==OptStrategy.none
%         descDir = -totalGrad;
%     end

    render(img, mesh, extra, approx, totalGrad, salmap, renderparams);
    
    %% check convergence and either subdivide or stop
    % if energy hasn't dropped significantly since window iterations ago, then energy is 'flat'
    if false
        error
    if i > windowSize && energy(i-windowSize) - energy(i) < demandedEnergyDensityDrop*totalArea 
        if subdivcount <= subdivmax
            subdivtimestart = tic;
            %% handle bad sliver triangles
            badTriInds = getTrisToCollapse(mesh,extra.perTriangleRGBError,img);
            if numel(badTriInds)~=0
                %{
                figure; set(gcf,'color','w');
                ax1 = subplot(1,1,1); image(img); hold all; axis equal; axis off; 
                patch('faces',preMesh.T,'vertices',preMesh.X,'facecolor','none','facealpha',0,'edgecolor','cyan','linewidth',1.5)
                patch('faces',T(badTriInds,:),'vertices',X,'facecolor','none','facealpha',1,'edgecolor','red','linewidth',1.5)

                figure; set(gcf,'color','w');
                ax2 = subplot(1,1,1); image(img); hold all; axis equal; axis off; 
                patch('faces',mesh.T,'vertices',mesh.X,'facecolor','none','facealpha',0,'edgecolor','cyan','linewidth',1.5)
                
                linkprop([ax1 ax2],{'xlim','ylim','zlim'});
                %}
                preMesh = mesh;
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
            if all(salmap(:)==1)
                subdivapprox = approx;
            end
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
            timedata.subdivtime(i) = toc(subdivtimestart);
            continue;
        else
            display('Optimization finished!');
            break;
        end
    end
    end
    
    %% find a good dt to use
    if dtstrat == DtStrategy.constrained
        dt = 1/max(vecnorm(descDir,2,2));
        error
        while any(getTriangleAreas(X+dt*descDir,T)<0)
            dt = dt / 2;
        end
    elseif dtstrat == DtStrategy.onepix
        
        dt = 5/max(vecnorm(gradToflatGrad(totalGrad),2,2));
    elseif dtstrat == DtStrategy.none
        if i==1; warning('constant dt is VERY not recommended.'); end
    end
    dts(i) = dt;

    %% update state
    X = X - dt * grad.Xgrad;
    edgeX = edgeX - dt * grad.edgeGrad;
    faceX = faceX - dt * grad.faceGrad;
    
    X = clipVerts(X,width,height);
    edgeX = permute(clipVerts(permute(edgeX,[1 3 2]),width,height),[1 3 2]);
    faceX = permute(clipVerts(permute(faceX,[1 3 2]),width,height),[1 3 2]);
end
    
render(img,mesh,colors,approx,[],[]);
title(sprintf('(iter:%d) (subdiviters:%d)', i, subdivcount-1));
    
if saveOut
    Xs{i} = X; Ts{i} = T;
    A = getframe(gcf);
    writeVideo(v, A.cdata);
    close(v);
    save([outpath 'XTs.mat'],'Xs','Ts')
    save([outpath 'timedata.mat'],'timedata')
end

figure; set(gcf,'color','w');
subplot(3,1,1); hold all; title('energy'); plot(energy(1:i-1)); 
for j=1:subdivmax; xline(subdiviters(j)); end
subplot(3,1,2); hold all; title('gradnorm'); plot(gradnorms(1:i-1)); 
subplot(3,1,3); hold all; title('dt'); plot(dts(1:i-1)); 

timedata.totalTime = toc(functionstarttime);

end

function flatGrad = gradToflatGrad(grad)
    flatGrad = [grad.Xgrad;...
        reshape(permute(grad.faceGrad,[1 3 2]),[],2);...
        reshape(permute(grad.edgeGrad,[1 3 2]),[],2),...
        ];
end
% 
% function [Xgrad, faceGrad, edgeGrad] = flatGradToGrad(flatGrad, mesh)
%     flatGrad(1:mesh.nX,:);
%     
%     flatGrad = [grad.Xgrad;...
%         reshape(permute(grad.faceGrad,[1 3 2]),[],2),...
%         reshape(permute(grad.edgeGrad,[1 3 2]),[],2),...
%         ];
% end
