function [extra, energy, grad] = polyComputeEnergy(img, mesh, integral1DNsamples, salmap, polyparams, colstrat)
    
    [extra, energy, colormaps] = polyComputeEnergyCore(img, mesh, integral1DNsamples, salmap, polyparams, colstrat);
    
    if nargout >= 3
        d = size(mesh.X,2);
        nT = mesh.nT;
        [~, prePerm] = XedgeXfaceXtoControlPoints(mesh);
        X0 = reshape(permute(prePerm.TX,[3 1 2]),[],2); % vertex positions per triangle. with redundancy.
        nFC = size(prePerm.FX,3); % number of face control points
        FaceX0 = prePerm.FX; % (nF d nFC)
        nEC = size(prePerm.EX,3)/3; % number of edge interior control points
        EdgeX0 = reshape(permute(reshape(prePerm.EX,nT,d,3,nEC),[3 1 2 4]),3*mesh.nT,d,nEC); % nnE d nEC. 
        newT = reshape(1:(3*nT),3,[])';
        newEdges = reshape(newT(:,[1 2 2 3 3 1])',2,[])';
        % soupMesh = PolyMeshFromXT(X0, newT, newEdges, EdgeX0, FaceX0, true, mesh.gdeg);
        
        %% get gradient w.r.t. corner vertices per triangle.
        eps = 5; % moves in unites of pixels
        flatXgrad = zeros(size(X0));
        perts = permute(reshape(eye(6)==1,3,2,6),[1 3 2]);
        for i=1:6
            pert = repmat(perts(:,i,:),1,nT,1);
            Xi = X0 + reshape(pert,[],2)*eps;
            meshI = PolyMeshFromXT(Xi, newT, newEdges, EdgeX0, FaceX0, true, mesh.gdeg);
            extraI = polyComputeEnergyCore(img, meshI, integral1DNsamples, salmap, polyparams, colstrat, colormaps);
            
            Xim = X0 - reshape(pert,[],2)*eps;
            meshIM = PolyMeshFromXT(Xim, newT, newEdges, EdgeX0, FaceX0, true, mesh.gdeg);
            extraIM = polyComputeEnergyCore(img, meshIM, integral1DNsamples, salmap, polyparams, colstrat, colormaps);
            fdiffI = (extraI.perTriangleRGBError - extraIM.perTriangleRGBError)/(2*eps);
            
            flatXgrad(pert) = fdiffI;
        end
        
        %% get gradient w.r.t. face vertices per triangle.
        flatFaceXgrad = zeros(size(FaceX0));
        perts = reshape(eye(nFC*d),d*nFC,d,nFC)==1;
        for i=1:d*nFC
            pert = perts(i,:,:);
            FaceXi = FaceX0 + pert*eps;
            meshI = PolyMeshFromXT(X0, newT, newEdges, EdgeX0, FaceXi, true, mesh.gdeg);
            extraI = polyComputeEnergyCore(img, meshI, integral1DNsamples, salmap, polyparams, colstrat, colormaps);
            
            FaceXim = FaceX0 - pert*eps;
            meshIM = PolyMeshFromXT(X0, newT, newEdges, EdgeX0, FaceXim, true, mesh.gdeg);
            extraIM = polyComputeEnergyCore(img, meshIM, integral1DNsamples, salmap, polyparams, colstrat, colormaps);
            
            fdiffI = (extraI.perTriangleRGBError - extraIM.perTriangleRGBError)/(2*eps);
            flatFaceXgrad(repmat(pert,nT,1,1)) = fdiffI;
        end
        
        %% get gradient w.r.t. edge vertices per triangle.
        flatEdgeXgrad = zeros(size(EdgeX0));
        rsEdgeX0 = reshape(EdgeX0,3,nT,2,nEC);
        perts = permute(reshape(eye(3*nEC*d),3*d*nEC,3,d,nEC),[2 1 3 4]) == 1;
        for i=1:3*d*nEC
            pert = perts(:,i,:,:);
            EdgeXi = reshape(rsEdgeX0 + pert*eps, 3*nT,d,nEC);
            meshI = PolyMeshFromXT(X0, newT, newEdges, EdgeXi, FaceX0, true, mesh.gdeg);
            extraI = polyComputeEnergyCore(img, meshI, integral1DNsamples, salmap, polyparams, colstrat, colormaps);
            
            EdgeXim = reshape(rsEdgeX0 - pert*eps, 3*nT,d,nEC);
            meshIM = PolyMeshFromXT(X0, newT, newEdges, EdgeXim, FaceX0, true, mesh.gdeg);
            extraIM = polyComputeEnergyCore(img, meshIM, integral1DNsamples, salmap, polyparams, colstrat, colormaps);
            
            fdiffI = (extraI.perTriangleRGBError - extraIM.perTriangleRGBError)/(2*eps);
            flatEdgeXgrad(repmat(pert,1,nT,1,1)) = fdiffI;
        end
        
        %% compile gradients together to form one aggregate gradient
        grad.faceGrad = flatFaceXgrad;
        
        grad.edgeGrad = zeros([mesh.nE*nEC, 2]);
        if nEC > 0
            [idxs, dists] = knnsearch(reshape(permute(EdgeX0,[1 3 2]),[],2), reshape(permute(mesh.edgeX,[1 3 2]),[],2),'K',2);
            assert(size(idxs,2)==2); % otherwise, it means there's fewer than 2 points to knn on. so geometry is low degree of mesh has less than 2 edges.
            interiorEdgeCPs = sum(abs(dists),2) < .01;
            tmp = reshape(permute(flatEdgeXgrad,[1 3 2]),[],2);
            grad.edgeGrad(interiorEdgeCPs,:) = tmp(idxs(interiorEdgeCPs,1),:) + tmp(idxs(interiorEdgeCPs,2),:); % add gradients from adjacent triangles
            grad.edgeGrad(~interiorEdgeCPs,:) = tmp(idxs(~interiorEdgeCPs,1),:); % boundary has only one edge contribution
        end
        grad.edgeGrad = permute(reshape(grad.edgeGrad,mesh.nE, nEC, 2),[1 3 2]); % nE, d, nEC
        
        flatXgrad = reshape(flatXgrad,3,nT,2);
        Tinds = reshape(mesh.T',[],1);
        grad.Xgrad = [accumarray(Tinds, reshape(flatXgrad(:,:,1),[],1), [mesh.nX 1]),...
                      accumarray(Tinds, reshape(flatXgrad(:,:,2),[],1), [mesh.nX 1])];
                  
      %% ensure boundary is slip conditions
      grad = polySlipConditions(mesh, grad);
    end
    
end