clear all; close all;

f = figure; hold all; set(gcf,'color','w');  
f.Position = [197.0000  399.4000  450.4000  299.2000];
title('mesh topology update');
xlabel('number of triangles'); ylabel('seconds')
x1s = []; y1s = [];
x2s = []; y2s = [];
for i=1:5
    load(sprintf('timedata%d.mat',i));
    keepinds = find(timedata.subdivtime~=0)';
    x1s = [x1s; timedata.nTs(keepinds)'];
    y1s = [y1s; timedata.energycomptime(keepinds)'];
end
for i=6:10
    load(sprintf('timedata%d.mat',i));
    keepinds = find(timedata.subdivtime~=0)';
    x2s = [x2s; timedata.nTs(keepinds)'];
    y2s = [y2s; timedata.energycomptime(keepinds)'];
end
scatter(x1s, y1s,'r.')
scatter(x2s, y2s,'g.')
legend({'linear','constant'})
exportgraphics(f,'topruntime.pdf')



f = figure; hold all; set(gcf,'color','w')
f.Position = [197.0000  399.4000  450.4000  299.2000];
title('mesh geometry update');
xlabel('number of triangles'); ylabel('seconds')
x1s = []; y1s = [];
x2s = []; y2s = [];
for i=1:5
    load(sprintf('timedata%d.mat',i));
    x1s = [x1s; timedata.nTs(:)];
    y1s = [y1s; timedata.energycomptime(:)];
end
for i=6:10
    load(sprintf('timedata%d.mat',i));
    x2s = [x2s; timedata.nTs(:)];
    y2s = [y2s; timedata.energycomptime(:)];
end
scatter(x1s, y1s,'r.')
scatter(x2s, y2s,'g.')
legend({'linear','constant'})
exportgraphics(f,'geomruntime.pdf')


