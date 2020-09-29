f1 = openfig('apple_energy_gradesc.fig')
f2 = openfig('apple_energy_rmsprop.fig')

x2 = f2.Children.Children(1).XData
y2 = f2.Children.Children(1).YData

x1 = f1.Children.Children(1).XData
y1 = f1.Children.Children(1).YData

f3 = figure; hold all; set(gcf,'color','w') 
f3.Position = [ 1.0000    1.0000  437.6000  195.2000];
xlim([0 228])
plot(x1,y1,'r');
plot(x2,y2,'g');
legend({'Gradient Descent','RMSProp'});
xlabel('Iterations');
ylabel('Energy');
ylim([min([y1(:);y2(:)]) max([y1(:);y2(:)])])
exportgraphics(f3,'jointEnergyFigs.pdf');