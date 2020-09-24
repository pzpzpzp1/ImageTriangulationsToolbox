function esMesh = buildTriSplitMesh(mesh)
    if nargin == 0
%         [X,T] = initialGridMesh(10, 10, 5, 1);
%         mesh = MeshFromXT(X,T);
        
        load 'temp2.mat';
    end
    X = mesh.X; T = mesh.T;
    
    newX = (X(mesh.edges(:,1),:)+X(mesh.edges(:,2),:))/2;
    Xnew = [X; newX];
    newvinds = (mesh.nX+1):size(Xnew,1);
    
    v1 = mesh.T(:,1);     v2 = mesh.T(:,2);    v3 = mesh.T(:,3);
    v4 = newvinds(mesh.triangles2edges(:,1))';
    v5 = newvinds(mesh.triangles2edges(:,2))';
    v6 = newvinds(mesh.triangles2edges(:,3))';
    
    Tnew = [v1 v4 v6; v4 v2 v5; v6 v5 v3; v6 v4 v5];
    Tnew = reshape(permute(reshape(Tnew,mesh.nT,4,3),[2 1 3]),[],3);
    
    flipped = find(getTriangleAreas(Xnew,Tnew)<=0);
    if numel(flipped)~=0
        newmin = min(getTriangleAreas(Xnew,Tnew));
        error(sprintf('Input mesh has triangle of minimum area %.16f. When subdivided this became %.16f. Error will occur.',min(mesh.triAreas), newmin));
    end
    
    esMesh = MeshFromXT(Xnew,Tnew);
    
    %{
    flipped = find(getTriangleAreas(Xnew,Tnew)<=0);
    preflipped = ceil(flipped/4);
    
    % before
    figure; 
    ax1=subplot(1,2,1); hold all;
    % axis equal; 
    patch('faces',T,'vertices',X,'facecolor','green','facealpha',.1)
    patch('faces',T(preflipped,:),'vertices',X,'facecolor','cyan','edgecolor','cyan')
    
    % after
        
    ax2=subplot(1,2,2); hold all;
    % axis equal; 
    patch('faces',Tnew,'vertices',Xnew,'facecolor','green','facealpha',.1)
    patch('faces',Tnew(flipped,:),'vertices',Xnew,'facecolor','cyan','edgecolor','cyan')
    
    linkprop([ax1 ax2],{'xlim','ylim','zlim'});
    %}
    
    
end