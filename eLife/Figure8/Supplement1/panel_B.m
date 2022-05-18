files = dir("./ant*");
skip = 0;  % number of files to skip (set to 2 if using the entire dir to remove "." and ".."
dirFlags = [files.isdir];
subFolders = files(dirFlags);
data = struct;
fignum = 1;
for i = (skip + 1):length(subFolders)
    [M0, M, V, U] = extract_data(convertCharsToStrings(subFolders(i).name));
    Uint = compute_U(M);
    t = 3*(1:190); %3*(1:size(Vint,1));
    data(i-skip).name = convertCharsToStrings(subFolders(i).name);
    data(i-skip).M = M;
    data(i-skip).V = V;
    data(i-skip).U = U;
    data(i-skip).t = 3*(1:size(Vint,1));
    %data(i-skip).Vint = Vint;

    figure(fignum);
    set(gcf,'color','w');
    ax = gca;
    plot(ax,t,Uint(1:190,1), LineWidth=1.75, Color='red');
    hold on
    plot(ax,t,-Uint(1:190,2), LineWidth=1.75, Color='#007f00');
    plot(ax,t,Uint(1:190,4), LineWidth=1.75, Color='blue');

    ax.Title.String = strcat('Domain : ', num2str(fignum), ' (', strrep(data(i-skip).name, '_',' '), ')');
    ax.Box = "on";
    ax.XLabel.String = 'time [min]';
    ax.XLabel.FontSize = 14;
    ax.YLabel.String = "$\hat{U}(t)$";
    ax.YLabel.Interpreter ='latex';
    ax.YLabel.FontSize = 16;
    %ax.YLim = [-0.45, 0.25];
    ax.Units = "normalized";
    ax.FontSize = 14; 
    ax.LineWidth = 1.0;
    legend(["$\varepsilon^{U}_{xx}$", "$\varepsilon^{U}_{xy}$", "$\varepsilon^{U}_{yy}$"], FontSize=18, Location="southwest", Interpreter='latex', Orientation='horizontal')
    legend boxoff 
    export_fig(sprintf('U_ant_domain_%d', fignum), '-pdf')
    fignum = fignum + 1;
end


