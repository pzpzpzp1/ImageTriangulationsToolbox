% colors is (3 vertices) (3 rgb) (nT triangles)
function [energy, colors, gradient] = linearComputeEnergy(img, mesh, integral1DNsamples)
    
    X = mesh.X; T = mesh.T; nT = size(T,1);
    % generate sample locations in barycentric coords
    ws = getBarycentricSamplingWeights(integral1DNsamples);
    samplePoints = getSamplePointsFromBarycentricWeights(ws, X, T);
    n = size(ws,1); 
   
    % perform sampling
    sampleVals = double(sampleImage(img, samplePoints)); % n, nT, 3
    
    % compute Lj: nT, 3, 3 
    % triangles, rgb colors, basis funcs
    Lj = squeeze(sum(sampleVals .* reshape(ws,n,1,1,3),1)).*mesh.triAreas/n;
    
    [K, Ki] = loadMassMatrix(1);
    Ki_T = repmat(reshape(Ki,1,3,3),nT,1,1);
    colors = multiprod(permute(Ki_T./(2*mesh.triAreas),[2 3 1]), permute(Lj,[3 2 1]));
    
    % compute actual energy.
    linearColor = permute(reshape(ws * reshape(colors,3,[]), n, 3, nT),[1 3 2]);
    energyPerColorChannel = mesh.triAreas'*squeeze(vecnorm(linearColor - sampleVals,2,1).^2)/n;
    energy = sum(energyPerColorChannel);
    
    % todo: optimize. should be slightly redundant with above.
    if nargout >= 3
        %% compute gradient
        X = mesh.X; T = mesh.T; edges = mesh.edges;
        nX = size(X,1); nT = size(T,1); nE = size(edges,1);

        % generate area samples of f
        ws = getBarycentricSamplingWeights(integral1DNsamples); n = size(ws,1); % number sample points
        samplePoints = getSamplePointsFromBarycentricWeights(ws, X, T); % n x nT x 2
        f_triangle = double(sampleImage(img, samplePoints)); % n x nT x 3

        % int_T_fphi_dA: (nT) (3 rgb) (3 phi_j)
        int_T_fphi_dA = squeeze(sum(f_triangle.*reshape(ws,n,1,1,3),1));

        % int_vn_dl: (nT) (6 = (2 xy) x (3 vertices))
        int_vn_dl = reshape(mesh.dAdt,[],6);

        % constantGradients: triangles, xy, phi
        constantGradients = mesh.dAdt./mesh.triAreas;
        % scratch: dphi_dt(:,:,phind,xyind,vertind) = constantGradients(:,xyind,phind).*ws(:,vertind);
        % dphi_dt: triangles, samples, phi, xy, verts
        dphi_dt = reshape(constantGradients,nT,1,1,2,3).*reshape(ws,1,n,3);

        % vn is (n, nT, 3, 6). 3 for number of edges. 6 for number of velocity values.
        vndl = sampleVdotN_dl(mesh, integral1DNsamples);

        % generate edge samples of f
        edgeSamplePoints = getEdgeSamplePoints(mesh,integral1DNsamples);
        f_tri_edges = sampleImage(img, edgeSamplePoints);

        % generate edge samples of phi
        % phi_edges_oneTri: (samples per edge) (3 edges) (3 phis)
        edgeWs = linspace(0,1,integral1DNsamples);
        phi_edges_oneTri = zeros(integral1DNsamples,3,3);
        phi_edges_oneTri(:,1,1) = fliplr(edgeWs);
        phi_edges_oneTri(:,2,2) = fliplr(edgeWs);
        phi_edges_oneTri(:,3,3) = fliplr(edgeWs);
        phi_edges_oneTri(:,3,1) = edgeWs;
        phi_edges_oneTri(:,1,2) = edgeWs;
        phi_edges_oneTri(:,2,3) = edgeWs;
        phi_edges = repmat(reshape(phi_edges_oneTri,1,integral1DNsamples,3,3),nT,1,1,1,1);

        % int_fvn_dl: (nT) (v motions) (rgb) (phi)
        int_fvn_dl = squeeze(sum(reshape(phi_edges,nT,integral1DNsamples,3,1,1,3).*permute(vndl,[2 1 3 4]).*reshape(double(f_tri_edges),nT,integral1DNsamples,3, 1, 3),[2 3]));

        % Build gradient preparation mat. Still vectorized per triangle. Will be re-indexed to lie on vertices.
        % todo: build gradprep
        [~, colors] = linearComputeEnergy(img, mesh, integral1DNsamples);
        [K, Ki] = loadMassMatrix(1);
        dcdt = Ki/2 * (dLdt/mesh.triAreas - Lj.*mesh.dAdt/mesh.triAreas.^2)
        % todo: build gradprep



        % accumulate per triangle gradient values to vertices.
        vertGrad = zeros(nX,6);
        for i=1:6
            vertGrad(:,i) = accumarray(mesh.T(:), reshape(gradPrep(:,:,i),[],1));
        end
        vertGrad = reshape(vertGrad,nX,2,3);

        % sum over rgb channels
        gradient = sum(vertGrad,3);

        % account for mesh boundary
        gradient = slipConditions(mesh, gradient);
    end
end