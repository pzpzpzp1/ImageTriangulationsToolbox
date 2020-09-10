function mesh = MeshFromXT(X,T)
    mesh.X = X;
    mesh.T = T;
    
    v1 = X(T(:,1),:);  v2 = X(T(:,2),:);  v3 = X(T(:,3),:);
    e12 = [v1-v2]; e12(1,3) = 0; e23 = [v2-v3]; e23(1,3) = 0;
    mesh.triAreas = (cross(e12,e23)/2)*[0 0 1]';
    
    
    
end