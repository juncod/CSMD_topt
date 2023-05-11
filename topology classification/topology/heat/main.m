clear; clc; close all;
tic
addpath('data');
[NODE, ELEM] = inp_('main128.inp');
BC_num = '32'; % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

penal = 3;
plotCut = 0.1;

name_num = 1;
for rand_num = 1:2
    for rmin_num = 1:2
        rmin = rmin_num * 0.4 - 0.2;
        [Hs, H] = prepare_filter(rmin, NODE, ELEM);
        for volfrac_num = 1:2
            volfrac = volfrac_num * 0.1 + 0.2;
            saveFileName = strcat('heat_',num2str(BC_num), '_', num2str(name_num));
            [x,BC_map] = topology(NODE, ELEM, volfrac, penal, saveFileName, plotCut, Hs, H,rand_num);
            name_num = name_num+1;
        end
    
    end
%    saveBCname = strcat('BC_','heat_train_',num2str(BC_num));
%    saveFunction_BCmap(saveBCname, BC_map)
end
disp('Finished-----------------')
toc
beep; close all;
%% Function
function [x,BC_map] = topology(NODE, ELEM, volfrac, penal, saveFileName, plotCut, Hs, H,rand_num)
    k0 = 1; kmin = 1e-3;
    nele = length(ELEM);
    Ue = zeros(4, 1);
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

    % f1 = figure; colormap(gray); axis equal; axis tight; axis off; caxis([-1 0]);
    f3 = figure; colormap(gray); axis equal; axis tight; axis off; caxis([-1 0]);
    % Start Iteration
    while change > 0.001 && iter < maxiter
        iter = iter + 1;
        if iter <= 15, gf = 1; else gf = min(1.5, 1.01 * gf); end
        xold = x;
        % FEA analysis
        [U, KE,BC_map] = FE_(NODE, ELEM, x, penal, k0, kmin, BC_mat);
        % Objective Function & Sensitivity Analysis
        c = 0;
        dc = zeros(1, nele);

        for i = 1:nele

            for j = 1:4
                Ue(j, 1) = U(ELEM(i, j));
            end

            c = c + (kmin + x(i) ^ penal * (k0 - kmin)) * Ue' * KE(:, :, i) * Ue;
            dc(i) = -penal * x(i) ^ (penal - 1) * (k0 - kmin) * Ue' * KE(:, :, i) * Ue;
        end

        % Filtering of Sensitivities
        dcn(:) = H * (x(:) .* dc(:)) ./ Hs ./ max(1e-3, x(:));
        % Design Update by the Optimality Criteria Method
        [x] = OC_(ELEM, x, volfrac, dcn, gf);
        % Print Results
        change = max(max(abs(x - xold)));
        disp([' It.: ' sprintf('%4i', iter) ' Obj.: ' sprintf('%10.4f', c) ...
                  ' Vol.: ' sprintf('%6.3f', sum(sum(x)) / (nele)) ...
                  ' ch.: ' sprintf('%6.3f', change)])
        % Plot Density
        plot_x = -ceil(max(0, x(:) - plotCut));
    end
    % figure(f1);
    % patch('Faces', ELEM, 'Vertices', NODE, 'FaceVertexCData', -x', 'FaceColor', 'flat', 'LineStyle', 'none');
    figure(f3);
    patch('Faces', ELEM, 'Vertices', NODE, 'FaceVertexCData', plot_x, 'FaceColor', 'flat', 'LineStyle', 'none');
    saveFunction(saveFileName, x, f3, plotCut)
end

function [U, KE,BC_map] = FE_(NODE, ELEM, x, penal, k0, kmin, BC_mat)
    n_node = length(NODE);
    nele = length(ELEM);
    dof = 1; sdof = n_node * dof;
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

    x_middle_up_section = abs(NODE(:, 1) - max(NODE(:, 1))/2) > 1e-6;
    x_middle_low_section = abs(NODE(:, 1) - max(NODE(:, 1))/2) < 1e-6;
    y_middle_low_section = abs(NODE(:, 2) - max(NODE(:, 1))/2) < 1e-6;

    BC = zeros(n_node, 1);
    BC_N = zeros(n_node, 1);
    BC(x_middle & y_lower, 1) = 1;
    BC(x_upper & y_middle, 1) = 1;
    %BC(x_middle_up & y_middle_up, 1) = 1; % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    BC_id = find(reshape(BC', [], 1));
    freedofs = find(reshape(~BC', [], 1));

    BC_N(x_lower, 1) = 1;  % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    BC_N(y_upper, 1) = 1;
%     BC_N(y_lower, 1) = 1;
%     BC_N(x_upper, 1) = 1;
    BC_Nid = find(reshape(BC_N', [], 1));
    %% Force condition %%
    F = zeros(sdof, 1);

%     F(:) = -0.025;
    F(BC_Nid(:)) = -1;  % => q, 열이 빠지는 곳

    [BC_nodes,~] = find(BC(:,1));
    [BC_N_nodes,~] = find(BC_N);
    for i = 1: size(BC_nodes,1)
        [elem_id,~] = find(ELEM == BC_nodes(i));
        BC_mat(elem_id) = -1;
    end
    for i = 1: size(BC_N_nodes,1)
        [elem_id,~] = find(ELEM == BC_N_nodes(i));
        BC_mat(elem_id) = 1;
    end
    [xlen,~] = find(NODE(:,1) == min(NODE(:,1)));
    [ylen,~] = find(NODE(:,2) == min(NODE(:,2)));
    BC_map = flip(reshape(BC_mat,size(xlen,1)-1,size(ylen,1)-1)');

    %% Solve Stiffness %%
    v = 0.3; h = 25;
    D = 1 / (1 - v ^ 2) * [1 v 0; v 1 0; 0 0 (1 - v) / 2];
    [KE, edofMat] = stiffness_(NODE, ELEM, D, h); % K = global stiffness, KE = local stiffness
    iK = reshape(kron(edofMat, ones(4, 1))', 16 * nele, 1);
    jK = reshape(kron(edofMat, ones(1, 4))', 16 * nele, 1);
    sK = zeros(16 * nele, 1);

    for i = 1:nele
        sK(16 * (i - 1) + 1:16 * i) = reshape(KE(:, :, i) * (kmin + x(i) ^ penal * (k0 - kmin)), 16, 1);
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
        E = intersect(ex, ey);

        for j = 1:length(E)
            k = k + 1;
            A = [m_elem(i, 1), m_elem(i, 2)];
            B = [m_elem(E(j), 1), m_elem(E(j), 2)];
            iH(k) = i;
            jH(k) = E(j);
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
