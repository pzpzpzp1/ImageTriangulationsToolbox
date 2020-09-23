% computes energy with and without tri split. tris that decrease in energy more when split will have higher score.
% score: [nT]. Higher scores for tris that should be split. 
function score = getTriSplitScore(mesh,img,approx,integral1DNsamples, salmap)

% build triangle soup of all possible divisions
esMesh = buildTriSplitMesh(mesh)
% compute energy per triangle of the soup
esExtra = approx.computeEnergy(img, esMesh, integral1DNsamples, salmap);
% compute the no split energy. this may be redundant with previous computation, but is quite fast and makes the function more standalone.
extra = approx.computeEnergy(img, mesh, integral1DNsamples, salmap);

% assemble the edge division comparisions
brokenTriEnergies = sum(permute(reshape(esExtra.perTriangleRGBError,4,[],3),[2 3 1]),3);
originalTriEnergies = extra.perTriangleRGBError;
score = sum(originalTriEnergies - brokenTriEnergies,2);

% error per triangle is computed including shared edges TWICE. 
% this means it's possible for the subdivision energy to be greater than the original. 
% this error goes away as integral samples increases, and you can clip any remaining outliers to 0.
score(score<0)=0;

end

%{
figure; image(img); hold all; axis equal; 
patch('vertices',mesh.X,'faces',mesh.edges(score>=0,:),'facecolor','none','edgecolor','green','linewidth',1);
patch('vertices',mesh.X,'faces',mesh.edges(score<0,:),'facecolor','none','edgecolor','k','linewidth',3);
patch('vertices',mesh.X,'faces',mesh.edges(score<0,:),'facecolor','none','edgecolor','w','linewidth',1);

figure;histogram(energyDrop(:))
%}