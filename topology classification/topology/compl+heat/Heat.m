function [c_heat, dc_heat] = Heat(x, NODE, ELEM, penal)

    nele = length(ELEM);
    n_node = length(NODE);
    dof = 1; sdof = n_node * dof;

    %%%%%%%%%%%%%%%%%%%% Boundary Condition %%%%%%%%%%%%%%%%%%%%
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

%     x_point1 = 2.5; y_point1 = 2.5;
%     node_point1 = abs(NODE(:, 1) - x_point1) < 1e-6 & abs(NODE(:, 2) - y_point1) < 1e-6;
%     x_point2 = 7.5; y_point2 = 7.5;
%     node_point2 = abs(NODE(:, 1) - x_point2) < 1e-6 & abs(NODE(:, 2) - y_point2) < 1e-6;
%     x_point3 = 6; y_point3 = 5;
%     node_point3 = abs(NODE(:, 1) - x_point3) < 1e-6 & abs(NODE(:, 2) - y_point3) < 1e-6;
%     x_point4 = 6; y_point4 = 4;
%     node_point4 = abs(NODE(:, 1) - x_point4) < 1e-6 & abs(NODE(:, 2) - y_point4) < 1e-6;
%     x_point5 = 5; y_point5 = 5;
%     node_point5 = abs(NODE(:, 1) - x_point5) < 1e-6 & abs(NODE(:, 2) - y_point5) < 1e-6;
    BC = zeros(n_node, 1);
    BC(x_upper & y_lower, :) = 1;
%     BC(x_middle_up & y_middle_up, :) = 1;
    % BC(node_point3 ,:) = 1;
    % BC(node_point4 ,:) = 1;
    % BC(node_point5 ,:) = 1;
    BCid = find(reshape(BC', [], 1));
    freedofs = find(reshape(~BC', [], 1));

    BC_N = zeros(n_node, 1);
    BC_N(x_lower, :) = 1;
    BC_Nid = find(reshape(BC_N', [], 1));
    %% Force condition %%
    F = zeros(sdof, 1); % => q, 열이 빠지는 곳
%     F(:) = -0.01;
    F(BC_Nid(:)) = -1;

    %%%%%%%%%%%%%%%%%%%% Solve Stiffness %%%%%%%%%%%%%%%%%%%%
    k0 = 1; kmin = 1e-3; v = 0.25; h = 25;
    D = 1 / (1 - v ^ 2) * [1 v 0; v 1 0; 0 0 (1 - v) / 2];
    [KE, edofMat] = stiffness_heat(NODE, ELEM, D, h); % K = global stiffness, KE = local stiffness
    iK = reshape(kron(edofMat, ones(4, 1))', 16 * nele, 1);
    jK = reshape(kron(edofMat, ones(1, 4))', 16 * nele, 1);
    sK = zeros(16 * nele, 1);

    for i = 1:nele
        sK(16 * (i - 1) + 1:16 * i) = reshape(KE(:, :, i) * (kmin + x(i) ^ penal * (k0 - kmin)), 16, 1);
    end

    K = sparse(iK, jK, sK); %K = (K+K')/2;
    U = zeros(sdof, 1); % displacement -> 2n-1 : u vector / 2n : v vector
    U(freedofs) = K(freedofs, freedofs) \ F(freedofs);

    %%%%%%%%%%%%%%%%%%%% Objective Function & Sensitivity Analysis %%%%%%%%%%%%%%%%%%%%
    Ue = zeros(4, 1);
    c_heat = 0;
    dc_heat = zeros(1, nele);

    for i = 1:nele

        for j = 1:4
            Ue(j, 1) = U(ELEM(i, j));
        end

        c_heat = c_heat + (kmin + x(i) ^ penal * (k0 - kmin)) * Ue' * KE(:, :, i) * Ue;
        dc_heat(i) = -penal * x(i) ^ (penal - 1) * (k0 - kmin) * Ue' * KE(:, :, i) * Ue;
    end
