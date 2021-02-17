
V = randn(6,2);
N = 100;

triMaps = poly_triangle_map(V);
[ws, interiorInds] = getBarycentricSamplingWeights(N); n = size(ws,1); 
[samplePoints, geometricGrad] = triMaps(ws); nT = 1;
geodets = reshape(get2DDeterminants(reshape(geometricGrad,[],2,2)),[],nT); % n nT
    
geodetsRS = geodets/max(abs(geodets));
xys = reshape(samplePoints,[],2);
ispos = geodets(:)>=0;
xypos = xys(ispos,:);
xyneg = xys(~ispos,:);
colpos = geodetsRS(ispos);
colneg = geodetsRS(~ispos);

%% visualize result
figurem; 
subplot(1,2,1); hold all; axis equal; colorbar;
xlim([min(xys(:,1)) max(xys(:,1))]);
ylim([min(xys(:,2)) max(xys(:,2))]);
scatter(xys(:,1),xys(:,2), 50, geodets(:),'filled')

subplot(1,2,2); hold all; axis equal; colorbar;
xlim([min(xys(:,1)) max(xys(:,1))]);
ylim([min(xys(:,2)) max(xys(:,2))]);
scatter(xypos(:,1),xypos(:,2), 10, 'g','filled')
scatter(xyneg(:,1),xyneg(:,2), 10, 'r','filled')



