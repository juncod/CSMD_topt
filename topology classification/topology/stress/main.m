clear; clc; close all;
addpath('FE'); addpath('MMA'); addpath('data');
[NODE, ELEM] = inp_('Job-L.inp');
BC_num = '14'; % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% pl:penal, q:stress relaxation, p:p-norm
pl = 3; q = 0.5; p = 10; plotCut = 0.25;

for rmin_num = 1:3
    rmin = rmin_num * 0.2;
    [Hs, H] = Prepare_filter(rmin, NODE, ELEM);

    for volfrac_num = 1:3
        volfrac = volfrac_num * 0.1 + 0.1;
        saveFileName = strcat(num2str(BC_num), '_', num2str(volfrac_num), '_', num2str(rmin_num));
        x = topology(NODE, ELEM, volfrac, pl, q, p, rmin, saveFileName, plotCut, Hs, H);
        close all;
    end

end

disp('Finished-----------------')
beep;

function x = topology(NODE, ELEM, volfrac, pl, q, p, rmin, saveFileName, plotCut, Hs, H)

    nele = length(ELEM);
    x = volfrac * ones(nele, 1);
    outeriter = 0;
    f1 = figure; colormap(gray); axis equal; axis tight; axis off; caxis([-1 0]);
    f2 = figure; colormap(jet); axis equal; axis tight; axis off; colorbar;
    f3 = figure; colormap(gray); axis equal; axis tight; axis off; caxis([-1 0]);
    % f5 = figure;
    % f6 = figure;

    % XY_elem = zeros(2,nele);
    % for i = 1:nele
    %     XY_elem(1,i) = sum(NODE(ELEM(i,:),1))/4;
    %     XY_elem(2,i) = sum(NODE(ELEM(i,:),2))/4;
    % end

    %%%%%%%%%%%%%%%%%% M M A Zone %%%%%%%%%%%%%%%%%%
    m = 1;
    epsimin = 0.0000001;
    n = length(x(:));
    xold1 = x;
    xold2 = x;
    xlb = 1e-3 * ones(n, 1);
    xub = 1 * ones(n, 1);
    xmin = xlb;
    xmax = xub;
    low = xlb;
    upp = xub;
    c = [1e4]';
    d = [0]';
    a0 = 1;
    a = [0]';
    raa0 = 0.0001;
    raa = 0.0001;
    raa0eps = 0.0000001;
    raaeps = 0.0000001;
    maxoutit = 120;
    kkttol = 0;
    nele = length(ELEM);
    x_his = zeros(nele, maxoutit);
    kktnorm = kkttol + 1;
    outit = 0;

    xIter = [1:maxoutit];
    volCons = zeros(1, maxoutit);
    strObjs = zeros(1, maxoutit);

    while outit < maxoutit
        outit = outit + 1;
        if outit <= 15, gf = 0.2; else gf = min(0.5, 1.01 * gf); end
        outeriter = outeriter + 1;
        [f0val, df0dx, fval, dfdx, MISES, sen, dens] = Stress_minimize(NODE, ELEM, x, Hs, H, pl, q, p, volfrac);
        %%%% The parameters low, upp, raa0 and raa are calculated:
        [low, upp, raa0, raa] = ...
        MMA_asymp(outeriter, n, x, xold1, xold2, xmin, xmax, low, upp, ...
            raa0, raa, raa0eps, raaeps, df0dx, dfdx);
        [xmma, ymma, zmma, lam, xsi, eta, mu, zet, s, f0app, fapp] = ...
            MMA_gcmmasub(m, n, outeriter, epsimin, x, xmin, xmax, low, upp, ...
            raa0, raa, f0val, df0dx, fval, dfdx, a0, a, c, d, gf);
        xold2 = xold1;
        xold1 = x;
        x = xmma;
        x_his(:, outit) = xmma;

        % PRINT RESULTS
        fprintf(' It.:%5i   P-norm Stress.:%11.4f   Vol.:%7.3f   MISES(max).:%11.4f   sens(max).:%11.4f   sens(min).:%11.4f    \n', outit, f0val, ...
        mean(dens(:)), max(MISES), max(df0dx), min(df0dx));

        % Plot Density
        plot_x = -ceil(max(0, dens(:) - plotCut));
        figure(f1);
        patch('Faces', ELEM, 'Vertices', NODE, 'FaceVertexCData', -dens, 'FaceColor', 'flat', 'LineStyle', 'none');
        figure(f2);
        patch('Faces', ELEM, 'Vertices', NODE, 'FaceVertexCData', MISES, 'FaceColor', 'flat', 'LineStyle', 'none'); caxis([min(MISES) max(MISES)]);
        figure(f3);
        patch('Faces', ELEM, 'Vertices', NODE, 'FaceVertexCData', plot_x, 'FaceColor', 'flat', 'LineStyle', 'none');

        %%%% The residual vector of the KKT conditions is calculated:
        % [residu,kktnorm,residumax] = ...
        % MMA_kktcheck(m,n,xmma,ymma,zmma,lam,xsi,eta,mu,zet,s, ...
        % xmin,xmax,df0dx,fval,dfdx,a0,a,c,d);
        %    outvector1 = [outeriter innerit x'];
        %    outvector2 = [f0val fval'];

        % volCons(1, outit) = mean(dens(:));
        % strObjs(1, outit) = max(MISES);

    end

    % figure(5);
    % plot(xIter, volCons(1, :));
    % xlabel('iter');
    % ylabel('volume fraction');

    % figure(6);
    % plot(xIter, strObjs(1, :));
    % xlabel('iter');
    % ylabel('Pnorm stress');

    saveFunction(saveFileName, x, f3, plotCut)

end
