clear; clc; close all;
tic
addpath('data');
[NODE, ELEM] = inp_('main128.inp');
BC_num = '59'; % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

penal = 3;
plotCut = 0.1;

name_num = 1;
for rand_num = 1:2
    for rmin_num = 1:2
        rmin = rmin_num * 0.4 - 0.2;
        [Hs, H] = prepare_filter(rmin, NODE, ELEM);
        for volfrac_num = 1:2
            volfrac = volfrac_num * 0.1 + 0.2;
            saveFileName = strcat('compl_',num2str(BC_num), '_', num2str(name_num));
            [x,BC_map] = topology(NODE, ELEM, volfrac, penal, saveFileName, plotCut, Hs, H,rand_num);
            name_num = name_num+1;
        end
    
    end
%     saveBCname = strcat('BC_','compl_test_',num2str(BC_num));
%     saveFunction_BCmap(saveBCname, BC_map)
end
disp('Finished-----------------')
toc
beep; close all;
%% Function
function [x,BC_map] = topology(NODE, ELEM, volfrac, penal, saveFileName, plotCut, Hs, H, rand_num)
    nele = length(ELEM);

    % XY_elem = zeros(2, nele);

    % for i = 1:nele
    %     XY_elem(1, i) = sum(NODE(ELEM(i, :), 1)) / 4;
    %     XY_elem(2, i) = sum(NODE(ELEM(i, :), 2)) / 4;
    % end

    % f4 = figure;

    Ue = zeros(8, 1);
    x(1:nele) = volfrac;
    BC_mat(1:nele) = 0;
    iter = 0;
    maxiter = 120;
    change = 1;
    if rand_num == 1
    else
        x = imnoise(x,'salt & pepper',1);
        maxiter = 150;
    end


%     f1 = figure; colormap(gray); axis equal; axis tight; axis off; caxis([-1 0]);
    f3 = figure; colormap(gray); axis equal; axis tight; axis off; caxis([-1 0]);
    
    c0_check = 0;
    c0 = 1;
    % Start Iteration
    while change > 0.001 && iter < maxiter
        iter = iter + 1;
        if iter <= 15, gf = 1; else gf = min(2, 1.01 * gf); end
        % RESTRICT DOMAIN
        % x(eZero_id) = 0.001;
        xold = x;
        % FEA analysis
        [U, KE,BC_map] = FE_(NODE, ELEM, x, penal, BC_mat);
        % Objective Function & Sensitivity Analysis
        c = 0;

        for i = 1:nele

            for j = 1:4
                Ue(2 * j - 1, 1) = U(2 * ELEM(i, j) - 1);
                Ue(2 * j, 1) = U(2 * ELEM(i, j));
            end

            c = c + (x(i) ^ penal) * Ue' * KE(:, :, i) * Ue;
            dc(i) = -penal * x(i) ^ (penal - 1) * Ue' * KE(:, :, i) * Ue;
        end

        % Filtering of Sensitivities
        dcn(:) = H * (x(:) .* dc(:)) ./ Hs ./ max(1e-3, x(:));

        % Design Update by the Optimality Criteria Method
        [x] = OC_(ELEM, x, volfrac, dcn, gf);
        % Print Results
        change = max(max(abs(x - xold)));
        
        if (isnan(c) || c0_check) == 0
            c0 = c;
            c0_check =1;
        end

        
        disp([' It.: ' sprintf('%4i', iter) ' Obj.: ' sprintf('%g', c/c0) ...
                  ' Vol.: ' sprintf('%6.3f', sum(sum(x)) / (nele)) ...
                  ' ch.: ' sprintf('%6.3f', change)])
        % Plot Density
        x(:) = real(x(:));
        plot_x = -ceil(max(0, x(:) - plotCut));

        % pause(1e-6);
        % figure(f4);
        % plot3(XY_elem(1, :), XY_elem(2, :), dcn, '.');
    end

%     figure(f1);
%     patch('Faces', ELEM, 'Vertices', NODE, 'FaceVertexCData', -x', 'FaceColor', 'flat', 'LineStyle', 'none');
    figure(f3);
    patch('Faces', ELEM, 'Vertices', NODE, 'FaceVertexCData', plot_x, 'FaceColor', 'flat', 'LineStyle', 'none');
    saveFunction(saveFileName, x, f3, plotCut)

end

function [U, KE,BC_map] = FE_(NODE, ELEM, x, penal, BC_mat)
    n_node = length(NODE);
    nele = length(ELEM);
    dof = 2; sdof = n_node * dof;
    %% Boundary Condition %%
    x_lower = abs(NODE(:, 1) - min(NODE(:, 1))) < 1e-6;
    x_upper = abs(NODE(:, 1) - max(NODE(:, 1))) < 1e-6;
    y_lower = abs(NODE(:, 2) - min(NODE(:, 2))) < 1e-6;
    y_upper = abs(NODE(:, 2) - max(NODE(:, 2))) < 1e-6;
    x_middle = abs(NODE(:, 1) - max(NODE(:, 1)) / 2) < 1e-6;
    y_middle = abs(NODE(:, 2) - max(NODE(:, 2)) / 2) < 1e-6; % % points of highest and lowest x,y % % points of highest and lowest x,y
    y_middle_low = abs(NODE(:, 2) - max(NODE(:, 2)) / 4) < 1e-6;
    y_middle_low_low = abs(NODE(:, 2) - max(NODE(:, 2)) / 8) < 1e-6;
    y_middle_low_up = abs(NODE(:, 2) - max(NODE(:, 2)) * 3/8) < 1e-6;
    x_middle_low = abs(NODE(:, 1) - max(NODE(:, 1)) / 4) < 1e-6;
    y_middle_up = abs(NODE(:, 2) - max(NODE(:, 2)) * 3/4) < 1e-6;
    x_middle_up = abs(NODE(:, 1) - max(NODE(:, 1)) * 3/4) < 1e-6;
    x_middle_up_low = abs(NODE(:, 1) - max(NODE(:, 1)) * 5/8) < 1e-6;
    x_middle_up_up = abs(NODE(:, 1) - max(NODE(:, 1)) * 7/8) < 1e-6;

    BC = zeros(n_node, 2);
    BC_N = zeros(n_node, 2);
%     BC(y_lower, :) = 1; % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    BC(y_upper, :) = 1;
    BC(x_lower, :) = 1;
%     BC(x_upper, :) = 1;
    BC_id = find(reshape(BC', [], 1));
    freedofs = find(reshape(~BC', [], 1));
    %% Force condition %%
    BC_N(x_middle & y_middle_low, 1) = 1; % % % % % % % % % % % % % % % % % % % % % % % % % % %
    BC_N(x_middle & y_middle_low, 2) = -1;
    BC_Nid = find(reshape(BC_N', [], 1));
    F = sparse(BC_Nid, 1, 1, sdof, 1);
    F(BC_Nid(:)) = 1;
    
    
    [BC_nodes,~] = find(BC(:,1));
    [BC_N_nodes,~] = find(BC_N);
    for i = 1: size(BC_nodes,1)
        [elem_id,~] = find(ELEM == BC_nodes(i));
        BC_mat(elem_id) = 1;
    end
    for i = 1: size(BC_N_nodes,1)
        [elem_id,~] = find(ELEM == BC_N_nodes(i));
        BC_mat(elem_id) = -1;
    end
    [xlen,~] = find(NODE(:,1) == min(NODE(:,1)));
    [ylen,~] = find(NODE(:,2) == min(NODE(:,2)));
    BC_map = flip(reshape(BC_mat,size(xlen,1)-1,size(ylen,1)-1)');

    %% Solve Stiffness %%
    E = 200e3; v = 0.3; h = 25;
    D = E / (1 - v ^ 2) * [1 v 0; v 1 0; 0 0 (1 - v) / 2];
    [KE, edofMat] = stiffness_(NODE, ELEM, D, h); % K = global stiffness, KE = local stiffness
    iK = reshape(kron(edofMat, ones(8, 1))', 64 * nele, 1);
    jK = reshape(kron(edofMat, ones(1, 8))', 64 * nele, 1);
    sK = zeros(64 * nele, 1);

    for i = 1:nele
        sK(64 * (i - 1) + 1:64 * i) = reshape(KE(:, :, i) * x(i)' .^ penal, 64, 1);
    end

    K = sparse(iK, jK, sK); %K = (K+K')/2;
    U = zeros(sdof, 1); % displacement -> 2n-1 : u vector / 2n : v vector
    U(freedofs) = K(freedofs, freedofs) \ F(freedofs);
end

%% mesh independcy filtering
function [Hs, H] = prepare_filter(rmin, NODE, ELEM)
    nele = length(ELEM);
    iH = zeros(nele * nele, 1); % 그냥 최대치 설정
    jH = zeros(size(iH));
    sH = zeros(size(iH));
    m_elem = zeros(nele, 2);

    for i = 1:nele
        m_elem(i, 1) = m_elem(i, 1) + sum(NODE(ELEM(i, :), 1)) / 4;
        m_elem(i, 2) = m_elem(i, 2) + sum(NODE(ELEM(i, :), 2)) / 4;
    end % => m_elem : element의 x좌표 y좌표

    k = 0;

    for i = 1:nele
        ex = find(abs(m_elem(:, 1) - m_elem(i, 1)) < rmin);
        ey = find(abs(m_elem(:, 2) - m_elem(i, 2)) < rmin);
        it = intersect(ex, ey);

        for j = 1:length(it)
            k = k + 1;
            A = [m_elem(i, 1), m_elem(i, 2)];
            B = [m_elem(it(j), 1), m_elem(it(j), 2)];
            iH(k) = i;
            jH(k) = it(j);
            sH(k) = max(0, rmin - norm(A - B));
        end

    end

    H = sparse(iH(1:k), jH(1:k), sH(1:k));
    Hs = sum(H, 2);
end

%% x new
function [xnew] = OC_(ELEM, x, volfrac, dcn, gf)
    l1 = 0; l2 = 1e5; move = 0.05;
    nele = length(ELEM);
    dv = ones(1, nele) / nele;

    while (l2 - l1 > 1e-6)
        lmid = 0.5 * (l1 + l2);
        xnew = max(0.001, max(x - move, min(1, min(x + move, (x .* sqrt(-dcn ./ dv ./ lmid)) .^ gf))));

        if sum(sum(xnew)) - volfrac * nele > 0
            l1 = lmid;
        else
            l2 = lmid;
        end

    end

end
