function badTriInds = getTrisToCollapse(mesh, perTriError, img)
    T = mesh.T; X = mesh.X; totalArea = size(img,1)*size(img,2);
    
    %% compute which triangles are bad
    a = mesh.triangleEdgeLengths(:,1);
    b = mesh.triangleEdgeLengths(:,2);
    c = mesh.triangleEdgeLengths(:,3);
    s = (a+b+c)/2;
    inRadii = sqrt(s.*(s-a).*(s-b).*(s-c))./s;

    triQual = inRadii./(max(mesh.triangleEdgeLengths')');
    areas = mesh.triAreas;
    triCosts = sqrt(sum(perTriError,2))./mesh.triAreas;

    % isBadTriShape = (triQual < prctile(triQual,10)) & (areas < prctile(areas,10)); stillbadtris = find(isBadTriShape);
    isBadTriShape = (triQual < .1) & (areas < .0009*totalArea); stillbadtris = find(isBadTriShape);
    triCosts = triCosts(isBadTriShape);
    badTriInds = stillbadtris(triCosts > prctile(triCosts,50));
    %{
    figure; 
    clf; image(img); hold all; axis equal; patch('faces',T(stillbadtris,:),'vertices',X,'facecolor','none','edgecolor','red')
    patch('faces',T(badTriInds,:),'vertices',X,'facecolor','none','edgecolor','blue')
    %}      
end        