function boost = getBoostMap(img, prevsal)
    f1 = figure('units','normalized'); image(img); 
    hold all;     axis equal;    
    set(gca,'XTickLabel',{},'YTickLAbel',{},'Box','on')
    if numel(prevsal)~=0
        prevsal = prevsal/max(prevsal(:));
        imh = image(prevsal*255);
        imh.AlphaData = prevsal*255/2;
    end
    
    
    BL = [0 0];
%     c1 = uicontrol('style','pushbutton','position',[BL(1) BL(2) 100 20],'string', 'Draw Polygon','callback',@drawpolygoncallback);
    c3 = uicontrol('style','pushbutton','position',[BL(1) BL(2)+20 100 20],'string', 'Draw Circle','callback',@drawcirclecallback);
%     polygons = {};
    circles = {};
%     function drawpolygoncallback(aa,bb) 
%         polygon = drawpolygon;
%         polygons{end+1} = polygon.Position;
%     end
    function drawcirclecallback(aa,bb) 
        circle = drawcircle;
        circles{end+1} = [circle.Center circle.Radius];
    end
    waitfor(f1);
    
    % img = img(1:2,1:3,:);
    boost = zeros(size(img,[1 2]));
    xvals = (1:size(img,2))-.5;
    yvals = (1:size(img,1))-.5;
    [X,Y] = meshgrid(xvals,yvals);
    % scatter(X(:),Y(:))
    
    for i=1:numel(circles)
        xyr = circles{i};
        
        ininds = (X(:)-xyr(1)).^2 + (Y(:)-xyr(2)).^2 - xyr(3)^2 < 0;
        boost(ininds) = boost(ininds) + 1;
    end
    
%     figure; 
%     image(img); hold all; axis equal;
%     imh = imagesc(boost); imh.AlphaData = .5;
    
end











