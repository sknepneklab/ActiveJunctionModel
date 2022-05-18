files = dir("./ant*");
skip = 0;  % number of files to skip (set to 2 if using the entire dir to remove "." and ".."
M0mode = "unit";   % check documenation for compute_Vint
dirFlags = [files.isdir];
subFolders = files(dirFlags);
data = struct;
Vmax = zeros(length(subFolders)-skip, 80);
Vmin = zeros(length(subFolders)-skip, 80);
Vmix = zeros(length(subFolders)-skip, 80);
fignum = 1;
for i = (skip + 1):length(subFolders)
    [M0, M, V, U] = extract_data(convertCharsToStrings(subFolders(i).name));
    Vint = compute_Vint(M0, M, V, M0mode);
    t = 3*(1:190); %3*(1:size(Vint,1));
    data(i-skip).name = convertCharsToStrings(subFolders(i).name);
    data(i-skip).M = M;
    data(i-skip).V = V;
    data(i-skip).U = U;
    data(i-skip).t = 3*(1:size(Vint,1));
    data(i-skip).Vint = Vint;

    figure(fignum);
    set(gcf,'color','w');
    ax = gca;
    plot(ax,t,3*Vint(1:190,1), LineWidth=1.75, Color='red');
    hold on
    plot(ax,t,-3*Vint(1:190,3), LineWidth=1.75, Color='#007f00');
    plot(ax,t,3*Vint(1:190,4), LineWidth=1.75, Color='blue');

    ax.Title.String = strcat('Domain : ', num2str(fignum), ' (', strrep(data(i-skip).name, '_',' '), ')');
    ax.Box = "on";
    ax.XLabel.String = 'time [min]';
    ax.XLabel.FontSize = 14;
    ax.YLabel.String = "$\int_0^t\hat{V}(t')dt'$";
    ax.YLabel.Interpreter ='latex';
    ax.YLabel.FontSize = 16;
    ax.YLim = [-0.45, 0.25];
    ax.Units = "normalized";
    ax.FontSize = 14; 
    ax.LineWidth = 1.0;
    legend(["$\varepsilon^{tot}_{xx}$", "$\varepsilon^{tot}_{xy}$", "$\varepsilon^{tot}_{yy}$"], FontSize=18, Location="southwest", Interpreter='latex', Orientation='horizontal')
    legend boxoff 
    export_fig(sprintf('ant_domain_%d', fignum), '-pdf')
    fignum = fignum + 1;
end


