clear; clc; close all;
addpath('FE'); addpath('data');
volfrac = 0.3; penal = 3; rmin = 0.3;

version = '1';
plotCut = 0.1;




tic
[NODE, ELEM] = inp_('main128.inp'); % % size: 10*10 (per 0.01)
[Hs, H] = Prepare_filter(rmin, NODE, ELEM);

compl_ratio = zeros(1, 21);
heat_ratio = zeros(1,21);
obj_ratio = linspace(0,1,21);
for i = 0:20
    alpha = i*(1/20);
    beta = 1 -alpha;
    ratio = sprintf('%1.2f',alpha);
    ratio_str = erase(ratio,'.');
    saveFileName = strcat(strcat(version,'_'),ratio_str);
    [c_ratio,h_ratio] = topology(NODE,ELEM, volfrac, penal, H, Hs, plotCut, saveFileName, alpha, beta);
    compl_ratio(i+1) = c_ratio;
    heat_ratio(i+1) = h_ratio;
    close all;
end
yyaxis left
xlabel('ObjectiveFunc ratio');
plot(obj_ratio,compl_ratio,'-o');
ylabel('Mechanical ratio');
yyaxis right
plot(obj_ratio,heat_ratio,'-o');
ylabel('Thermal ratio');
beep
toc

function [c_ratio,h_ratio] = topology(NODE,ELEM, volfrac, penal, H, Hs, plotCut, saveFileName, alpha, beta)
nele = length(ELEM);
x(1:nele) = volfrac;
iter = 0;
maxiter = 120;
change = 1;
f1 = figure; colormap(gray); axis equal; axis tight; axis off; caxis([-1 0]);
f3 = figure; colormap(gray); axis equal; axis tight; axis off; caxis([-1 0]);

x_graph = [1:maxiter];
y_compl = zeros(maxiter, 1);
y_heat = zeros(maxiter, 1);
y_complSens = zeros(maxiter, 1);
y_heatSens = zeros(maxiter, 1);


[c_compl, ~] = Compliance(x, NODE, ELEM, penal);
[c_heat, ~] = Heat(x, NODE, ELEM, penal);
gamma = c_compl/c_heat;
if alpha == 0
    gamma =1;
end
if beta == 0;
    gamma = 1;
end

while change > 0.001 && iter < maxiter
    iter = iter + 1;
    if iter <= 15, gf = 1; else gf = min(1.5, 1.01 * gf); end
    xold = x;
    [c_compl, dc_compl] = Compliance(x, NODE, ELEM, penal);
    [c_heat, dc_heat] = Heat(x, NODE, ELEM, penal);
    gamma =1;
    c_heat = gamma * c_heat;
    dc_heat = gamma * dc_heat;
    c = alpha * c_compl + beta * c_heat;
    dc = alpha * dc_compl + beta * dc_heat;
    % Filtering of Sensitivities
    dc(:) = H * (x(:) .* dc(:)) ./ Hs ./ max(1e-3, x(:));
    % Design Update by the Optimality Criteria Method
    [x] = OC_(ELEM, x, volfrac, dc, gf);
    % Print Results
    change = max(max(abs(x - xold)));
    y_compl(iter) = c_compl; y_heat(iter) = c_heat;
    y_complSens(iter) = max(abs(dc_compl)); y_heatSens(iter) = max(abs(dc_heat));
%     disp([' It.: ' sprintf('%4i', iter) ' Compl_ratio.: ' sprintf('%6.3f', c_compl/y_compl(1)) ' Heat_ratio.: ' sprintf('%6.3f', c_heat/y_heat(1)) ...
%               ' Vol.: ' sprintf('%6.3f', sum(sum(x)) / (nele)) ...
%               ' ch.: ' sprintf('%6.3f', change)])
%     disp([' Obj_compl.: ' sprintf('%10.4f', c_compl) ' Obj_heat.: ' sprintf('%10.4f', c_heat) ...
%               ' Sens_compl.: ' sprintf('%10.4f', max(abs(dc_compl))) ...
%               ' Sens_heat.: ' sprintf('%10.4f', max(abs(dc_heat)))])

    % Plot Density
    plot_x = -ceil(max(0, x(:) - plotCut));
    figure(f1);
    patch('Faces', ELEM, 'Vertices', NODE, 'FaceVertexCData', -x', 'FaceColor', 'flat', 'LineStyle', 'none');
    figure(f3);
    patch('Faces', ELEM, 'Vertices', NODE, 'FaceVertexCData', plot_x, 'FaceColor', 'flat', 'LineStyle', 'none');
    pause(1e-6);
end

    disp([' Alpha.: ' sprintf('%2.2f', alpha) ' Compl_ratio.: ' sprintf('%6.3f', c_compl/y_compl(1)) ' Heat_ratio.: ' sprintf('%6.3f', c_heat/y_heat(1)) ])

%     figure;
%     tiledlayout(2, 1)
%     ax1 = nexttile;
%     plot(ax1, x_graph, y_compl, x_graph, y_heat)
%     ax2 = nexttile;
%     plot(ax2, x_graph, y_complSens, x_graph, y_heatSens)
    saveFunction(saveFileName, x, f1, f3, plotCut)
    c_ratio = c_compl/y_compl(1);
    h_ratio = c_heat/y_heat(1);
end



