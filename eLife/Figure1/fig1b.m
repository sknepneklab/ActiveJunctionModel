Tst = 0.3;
k0 = 2/Tst;
FC = @(Text, beta) beta.*(1./(1+exp(k0*(Tst-Text))) - 0.5) - Text;

[Text,beta] = meshgrid(linspace(0,1,100), linspace(0,5,100));
%colormap(redblue)
set(gcf,'units','points','position',[0,0,450,300])
set(0, 'DefaultTextInterpreter', 'latex')
colormap(jet)
imagesc(Text(1,:), beta(:,1), FC(Text, beta));
cb = colorbar('XTick',-1.5:0.5:1.5);
cb.TickLabelInterpreter = 'latex';
%surf(Text, beta, FC(Text, beta));
caxis([-1.5,1.5]);
hold on
shading interp 
view(2);
[c,h] = contour(Text, beta, FC(Text, beta), [-1.4:0.2:-0.2 0.2:0.2:1.4], 'LineWidth', 0.75, 'LineColor', [1 1 1]);
hold on
clabel(c,h,'Color',[1 1 1],'Fontsize',14,'LabelSpacing',300, 'Interpreter','latex')
%clabel(c,h,'manual', 'Interpreter','latex')
contour(Text, beta, FC(Text, beta), [0 0], 'LineColor','k', 'LineWidth', 1)
%set(gcf,'units','points','position',[0,0,600,600])
ax = gca;
ax.YDir = 'normal';
ax.TickLabelInterpreter = 'latex';
ax.LineWidth = 0.5;
ax.FontSize = 18;
ax.XLabel.String = '$T_{ext}$';
ax.YLabel.String = '$\beta$';
ax.YTick = 0:1:5;
ax.XTick = 0:0.2:1;
print -painters -dpdf -r600 fig3a.pdf