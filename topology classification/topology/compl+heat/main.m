clear; clc; close all;
addpath('FE'); addpath('MMA'); addpath('data');
volfrac = 0.3; penal = 3; rmin = 0.3;

saveFileName = '0';
plotCut = 0.1;

[NODE,ELEM] = inp_('main.inp'); %% size: 10*10 (per 0.01)
[Hs,H]=Prepare_filter(rmin,NODE,ELEM);
nele = length(ELEM);
x(1:nele) = volfrac;
iter = 0;
maxiter = 30;
change = 1;
f1 = figure; colormap(gray); axis equal; axis tight; axis off; caxis([-1 0]);
f3 = figure; colormap(gray); axis equal; axis tight; axis off; caxis([-1 0]);

x_graph = [1:maxiter];
y_compl = zeros(maxiter,1);
y_heat = zeros(maxiter,1);
y_complSens = zeros(maxiter,1);
y_heatSens = zeros(maxiter,1);

while change > 0.001 && iter < maxiter
    iter = iter + 1;
    if iter <= 15, gf = 1; else gf = min(1.5,1.01*gf); end
    xold = x;
    [c_compl,dc_compl] = Compliance(x,NODE,ELEM,penal);
    [c_heat,dc_heat] = Heat(x,NODE,ELEM,penal);
    c = c_compl + c_heat;
    dc = dc_compl + dc_heat;
    % Filtering of Sensitivities
    dc(:) = H*(x(:).*dc(:))./Hs./max(1e-3,x(:));
    % Design Update by the Optimality Criteria Method
    [x] = OC_(ELEM,x,volfrac,dc,gf);
    % Print Results
    change = max(max(abs(x-xold)));
    disp([' It.: ' sprintf('%4i',iter) ' Obj.: ' sprintf('%10.4f',c) ...
   ' Vol.: ' sprintf('%6.3f',sum(sum(x))/(nele)) ...
    ' ch.: ' sprintf('%6.3f',change)])
    disp([' Obj_compl.: ' sprintf('%10.4f',c_compl) ' Obj_heat.: ' sprintf('%10.4f',c_heat) ...
   ' Sens_compl.: ' sprintf('%10.4f',max(abs(dc_compl))) ...
    ' Sens_heat.: ' sprintf('%10.4f',max(abs(dc_heat)))])
    % Plot Density
    plot_x = -ceil(max(0,x(:)-plotCut));
    figure(f1); 
    patch('Faces',ELEM,'Vertices',NODE,'FaceVertexCData',-x','FaceColor','flat','LineStyle','none'); 
    figure(f3); 
    patch('Faces',ELEM,'Vertices',NODE,'FaceVertexCData',plot_x,'FaceColor','flat','LineStyle','none'); 
    pause(1e-6);       
    y_compl(iter) = c_compl;   y_heat(iter) = c_heat;
    y_complSens(iter) = max(abs(dc_compl));   y_heatSens(iter) = max(abs(dc_heat));
end
tiledlayout(2,1)
ax1 = nexttile;
plot(ax1,x_graph,y_compl,x_graph,y_heat)
ax2 = nexttile;
plot(ax2,x_graph,y_complSens,x_graph,y_heatSens)
saveFunction(saveFileName,x,f1,f3,plotCut)
