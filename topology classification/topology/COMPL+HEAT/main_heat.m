clear; clc; close all;
addpath('FE'); addpath('data');
tic
[NODE, ELEM] = inp_('main128.inp');
BC_num = '888'; % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

penal = 3;
plotCut = 0.1;

name_num = 1;
% for rand_num = 1:2
    for rmin_num = 1:2
        rmin = rmin_num * 0.4 - 0.2;
        [Hs, H] = Prepare_filter(rmin, NODE, ELEM);
        for volfrac_num = 1:2
            volfrac = volfrac_num * 0.1 + 0.2;
            saveFileName = strcat('heat_test_',num2str(BC_num), '_', num2str(name_num));
            [x] = topology(NODE, ELEM, volfrac, penal, saveFileName, plotCut, Hs, H);
            name_num = name_num+1;
        end
    
    end
%     saveBCname = strcat('BC_','compl_test_',num2str(BC_num));
%     saveFunction_BCmap(saveBCname, BC_map)
% end
disp('Finished-----------------')
toc
beep; close all;

function [x] = topology(NODE, ELEM, volfrac, penal, saveFileName, plotCut, Hs, H)
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

while change > 0.001 && iter < maxiter
    iter = iter + 1;
    if iter <= 15, gf = 1; else gf = min(1.5, 1.01 * gf); end
    xold = x;
    [c_compl, dc_compl] = Compliance(x, NODE, ELEM, penal);
    [c_heat, dc_heat] = Heat(x, NODE, ELEM, penal);
%     c_heat = gamma * c_heat;
%     dc_heat = gamma * dc_heat;
%     c = alpha * c_compl + beta * c_heat;
%     dc = alpha * dc_compl + beta * dc_heat;
    c = c_heat; dc = dc_heat;
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

    disp([ ' Compl_ratio.: ' sprintf('%6.3f', c_compl/y_compl(1)) ' Heat_ratio.: ' sprintf('%6.3f', c_heat/y_heat(1)) ])

%     figure;
%     tiledlayout(2, 1)
%     ax1 = nexttile;
%     plot(ax1, x_graph, y_compl, x_graph, y_heat)
%     ax2 = nexttile;
%     plot(ax2, x_graph, y_complSens, x_graph, y_heatSens)

%     saveFunction(saveFileName, x, f1, f3, plotCut)
    c_ratio = c_compl/y_compl(1);
    h_ratio = c_heat/y_heat(1);
end
