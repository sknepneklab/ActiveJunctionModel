beta = 2.5;
Text = 0.5;
taum = 10; 
tauv = 10;
tspan = 0:0.01:200;

B = [0, 0.1, 0.2 0.228478681684888, 0.3];
tout = {};
yout = {};

for i = 1:numel(B) 
    options = odeset('MaxStep', 0.01);
    [t,y] = ode45(@(t,y) ajm_1d(t,y,beta,Text,B(i),taum,tauv), tspan, [0, 0.6, 1], options);
    tout{i} = t;
    yout{i} = y;
end

u = {};
m = {};
l0 = {};
l = {};
for i = 1:numel(B)
   u{i} = yout{i}(:,1);
   m{i} = yout{i}(:,2);
   l0{i} = yout{i}(:,3); 
   l{i} = u{i} + l0{i};
end

close all

cmap = colormap(jet);
for i = 1:numel(B)
    di = floor(256/numel(B));
    plot(tout{i},l{i},'-','LineWidth',1,'Color', cmap(i*di,:))
    hold on
end
ax = gca;
ax.TickLabelInterpreter = 'latex';
ax.LineWidth = 0.5;
ax.FontSize = 22;
ax.YTick = -0.2:0.2:1.2;
ax.XLabel.String = '$t/t^*$';
ax.YLabel.String = '$l/a$';
ax.XLim = [0,200];
ax.YLim = [-0.2,1.3];
figure
for i = 1:numel(B)
    di = floor(256/numel(B));
    plot(tout{i},m{i},'-','LineWidth',0.5,'Color', cmap(i*di,:))
    hold on
end
legend({'$0.0$','$0.1$', '$0.2$', '$FC\approx0.2285$', '$0.3$'},'Interpreter','latex','Location', 'south') 
ax = gca;
ax.TickLabelInterpreter = 'latex';
ax.LineWidth = 0.5;
ax.FontSize = 32;
ax.XLabel.String = '$t/t^*$';
ax.YLabel.String = '$m$';
ax.XLim = [0,200];
ax.YLim = [0.58,0.84];
ax.XTick = 0:50:200;

figure
for i = 1:numel(B)
    di = floor(256/numel(B));
    plot(m{i},u{i},'-','LineWidth',1,'Color', cmap(i*di,:))
    hold on
end

legend({'$0.0$','$0.1$', '$0.2$', '$FC\approx0.2285$', '$0.3$'},'Interpreter','latex','Location', 'southwest') 
ax = gca;
ax.TickLabelInterpreter = 'latex';
ax.LineWidth = 1;
ax.FontSize = 32;
ax.XLabel.String = '$m$';
ax.YLabel.String = '$u$';
ax.XLim = [0.56,0.84];
ax.YLim = [-0.35,0.35];
ax.XTick = 0.6:0.05:0.8;
ax.YTick = -0.3:0.1:0.3;