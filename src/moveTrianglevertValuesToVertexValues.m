% mesh needs mesh.T, mesh.X
% trivertValues: [nT 3 (m,n,...)]
% vertGrad: [nX (m,n,...)]
function vertVals = moveTrianglevertValuesToVertexValues(mesh, trivertValues)
    nX = size(mesh.X,1);
    nT = size(mesh.T,1);
    
    maxdim = ndims(trivertValues);
    mn = size(trivertValues,[3:maxdim]);
    
    flattened = reshape(trivertValues,nT,3,[]); 
    vertVals = zeros(nX,size(flattened,3));
    for i=1:size(flattened,3) % for loop meaning (m,n,...) should not be very big.
        vertVals(:,i) = accumarray(mesh.T(:), reshape(flattened(:,:,i),[],1));
    end
    vertVals = reshape(vertVals,[nX, mn]);
end