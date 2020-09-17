% Computes colors on vertices per triangle in two different ways. 
% First using all sampled colors over the triangle. 
% Second using a subset of sampled colors dicated by interiorInds.
% The second may be better for rendering purposes. Less sensitive to edge artifacts.
% mesh must have X and T
% samplePoints msut be   [n nT 2]
% f_triangle             [nT n 3]
% interiorInds           [n]
function [colors, colorsalt] = linearGetColorsFromSamples(mesh, samplePoints, f_triangle, interiorInds)
    X = mesh.X; T = mesh.T; nT = size(T,1); n = size(samplePoints,1);
    TX = reshape(X(T(:),:),nT,3,2);
    TX1 = TX; TX1(:,:,3)=1;
    
    %% least squares color using full sampling
    % LHS
    A1 = zeros(n,nT,3);    
    A1(:,:,1:2) = samplePoints; A1(:,:,3)=1;
    AtA = squeeze(sum(A1.*reshape(A1,n,nT,1,3),1));
    % RHS
    rhsR = squeeze(sum(f_triangle(:,:,1) .* permute(A1,[2 1 3]),2));
    rhsG = squeeze(sum(f_triangle(:,:,2) .* permute(A1,[2 1 3]),2));
    rhsB = squeeze(sum(f_triangle(:,:,3) .* permute(A1,[2 1 3]),2));
    % solve
    AtAi = multinv(permute(AtA,[2 3 1]));
    Rcoeffs = permute(sum(AtAi .* reshape(rhsR',1,3,nT),2),[3 2 1]);
    Gcoeffs = permute(sum(AtAi .* reshape(rhsG',1,3,nT),2),[3 2 1]);
    Bcoeffs = permute(sum(AtAi .* reshape(rhsB',1,3,nT),2),[3 2 1]);
    colors(:,:,1) = sum(TX1.*Rcoeffs,3);
    colors(:,:,2) = sum(TX1.*Gcoeffs,3);
    colors(:,:,3) = sum(TX1.*Bcoeffs,3);
    
    %% least squares color using limited sampling
    % LHS
    nm = sum(interiorInds);
    A1 = zeros(nm,nT,3);    
    A1(:,:,1:2) = samplePoints(interiorInds,:,:); A1(:,:,3)=1;
    AtA = squeeze(sum(A1.*reshape(A1,nm,nT,1,3),1));
    % RHS
    rhsR = squeeze(sum(f_triangle(:,interiorInds,1) .* permute(A1,[2 1 3]),2));
    rhsG = squeeze(sum(f_triangle(:,interiorInds,2) .* permute(A1,[2 1 3]),2));
    rhsB = squeeze(sum(f_triangle(:,interiorInds,3) .* permute(A1,[2 1 3]),2));
    % solve
    AtAi = multinv(permute(AtA,[2 3 1]));
    Rcoeffs = permute(sum(AtAi .* reshape(rhsR',1,3,nT),2),[3 2 1]);
    Gcoeffs = permute(sum(AtAi .* reshape(rhsG',1,3,nT),2),[3 2 1]);
    Bcoeffs = permute(sum(AtAi .* reshape(rhsB',1,3,nT),2),[3 2 1]);
    colorsalt(:,:,1) = sum(TX1.*Rcoeffs,3);
    colorsalt(:,:,2) = sum(TX1.*Gcoeffs,3);
    colorsalt(:,:,3) = sum(TX1.*Bcoeffs,3);
    
end