function scr = shortScatter(points, value, cutoff)

if nargin < 2
    value = ones(size(points,1),1);
end

if nargin < 4
    cutoff = mean(value) + std(value);
end

idx = find(value >= cutoff);
if numel(idx)==1
    return;
end


colors = value(idx);
colors = colors - min(colors) + 1e-6;
colors = colors / max(colors);

if size(points,2)==2
    scr = scatter(points(idx,1),points(idx,2), colors*10, colors, 'filled');
else
    scr = scatter3(points(idx,1),points(idx,2),points(idx,3), colors*10, colors, 'filled');
end

end