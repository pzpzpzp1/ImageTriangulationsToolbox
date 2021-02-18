function [energy, grad] = polyGetAreaEnergy(mesh, integral1DNsamples, salmap, fdiffMeshes)
    d=2;
    beta = @(x) 1./x;
    factor = beta(salmap);
    energy = sum(polyGetAreaEnergyCore(mesh, integral1DNsamples, factor));
    
    if nargout > 1
        assert(exist('fdiffMeshes','var')==1); % should be cached from previous fdiff energy
        soupmesh = fdiffMeshes.VertPlus(1);
        nT = mesh.nT;
        
        % vertex fdiffs
        flatXgrad = zeros(size(soupmesh.X));
        for i=1:6
            E1 = polyGetAreaEnergyCore(fdiffMeshes.VertPlus(i), integral1DNsamples, factor);
            E2 = polyGetAreaEnergyCore(fdiffMeshes.VertMinus(i), integral1DNsamples, factor);
            flatXgrad(repmat(fdiffMeshes.vertperts(:,i,:),1,soupmesh.nT,1)) = (E1-E2)/(2*fdiffMeshes.eps);
        end
        
        % face fdiffs
        flatFaceXgrad = zeros(size(soupmesh.faceX));
        nFC = size(flatFaceXgrad,3);
        for i=1:d*nFC
            E1 = polyGetAreaEnergyCore(fdiffMeshes.FacePlus(i), integral1DNsamples, factor);
            E2 = polyGetAreaEnergyCore(fdiffMeshes.FaceMinus(i), integral1DNsamples, factor);
            flatFaceXgrad(repmat(fdiffMeshes.faceperts(i,:,:),nT,1,1)) = (E1-E2)/(2*fdiffMeshes.eps);
        end
        
        % edge fdiffs
        flatEdgeXgrad = zeros(size(soupmesh.edgeX));
        nEC = size(flatEdgeXgrad,3);
        for i=1:3*d*nEC
            E1 = polyGetAreaEnergyCore(fdiffMeshes.EdgePlus(i), integral1DNsamples, factor);
            E2 = polyGetAreaEnergyCore(fdiffMeshes.EdgeMinus(i), integral1DNsamples, factor);
            flatEdgeXgrad(repmat(fdiffMeshes.edgeperts(:,i,:,:),1,nT,1,1)) = (E1-E2)/(2*fdiffMeshes.eps);
        end
        
        % re-aggegrate 
        grad.faceGrad = flatFaceXgrad;
        
        grad.edgeGrad = zeros([mesh.nE*nEC, 2]);
        if nEC > 0
            [idxs, dists] = knnsearch(reshape(permute(soupmesh.edgeX,[1 3 2]),[],2), reshape(permute(mesh.edgeX,[1 3 2]),[],2),'K',2);
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
                  
        grad = polySlipConditions(mesh, grad);
    end
end