
function [cp, prePerm] = XedgeXfaceXtoControlPoints(mesh, override)
    gdeg = mesh.gdeg;
    TX = permute(reshape(mesh.X(mesh.T',:),3,[],2),[2 3 1]);
    
    % get oriented edge control points per triangle
    nEC = max(gdeg-1,0); % edge control points.
    E1 = mesh.edgeX(mesh.triangles2edges(:,1),:,:);
    E2 = mesh.edgeX(mesh.triangles2edges(:,2),:,:);
    E3 = mesh.edgeX(mesh.triangles2edges(:,3),:,:);
    E123 = reshape(permute(cat(4,E1,E2,E3),[1 4 2 3]), mesh.nT*3,2,gdeg-1);
    E123(~mesh.triangles2edges_isForward,:,:) = flip(E123(~mesh.triangles2edges_isForward,:,:),3);
    EX = reshape(permute(reshape(E123,mesh.nT,3,2,nEC),[1 3 4 2]),mesh.nT,2,[]); % nT 2 3 nEC
    
    cp = cat(3,TX, EX, mesh.faceX);
    perm = 1:size(cp,3);
    if gdeg == 1
        perm = [3 2 1];
    elseif gdeg == 2
        perm = [3     5     2     6     4     1];
    elseif gdeg == 3
        perm = [3     6     7     2     9    10     4     8     5     1];
    elseif gdeg == 4
        perm = [3     7     8     9     2    12    13    14     4    11    15     5    10     6     1];
    elseif gdeg == 5
        perm = [3     8     9    10    11     2    15    16    17    18     4    14    19    20     5    13    21     6    12     7     1];
    else
        if ~(exist('override','var') && override)
            error('no permutation saved. obtain using getPerm manually, and hardcode it here. Though do you really need such high degree triangles?');
        end
    end
    prePerm.TX = TX;
    prePerm.FX = mesh.faceX;
    prePerm.EX = EX;
    cp = cp(:,:,perm);
end

% obtains permutation to take output of XedgeXfaceXtoControlPoints and rearrange it to match that needed for poly_triangle_map.
function perm = getPerm(mesh)
    gdeg = mesh.gdeg;
    triangleControlPoints = XedgeXfaceXtoControlPoints(mesh, 1);
    
    ws = getBarycentricSamplingWeights(gdeg+1); ws = reshape(ws,1,1,[],3);
    triangleControlPoints2 = mesh.X(mesh.T(:,1),:).*ws(1,1,:,1) + mesh.X(mesh.T(:,2),:).*ws(1,1,:,2) + mesh.X(mesh.T(:,3),:).*ws(1,1,:,3);
    ti = randi(mesh.nT); x1 = squeeze(triangleControlPoints(ti,:,:))'; x2 = squeeze(triangleControlPoints2(ti,:,:))';
    [perm,ds] = knnsearch(x1, x2); assert(norm(ds)<.00001); 
    perm = perm';
end