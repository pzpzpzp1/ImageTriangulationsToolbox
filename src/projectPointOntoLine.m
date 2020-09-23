% projects v1 onto the line specified by v2 -> v3
function v4 = projectPointOntoLine(v1,v2,v3)
    e12 = v1-v2;
    e23 = v3-v2;
    
    v4 = v2 + dot(e12,e23,2)./dot(e23,e23,2).*e23;
    
    %{
    i = 3; 
    figure; hold all; axis equal;
    scatter(v1(i,1), v1(i,2),'r','filled');
    v23 = [v2(i,:);v3(i,:)];
    plot(v23(:,1), v23(:,2),'g');
    scatter(v4(i,1), v4(i,2),'k','filled');
    %}
end