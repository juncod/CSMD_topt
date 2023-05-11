function [c_compl,dc_compl] = Compliance(x,NODE,ELEM,penal)

    nele = length(ELEM);
    n_node = length(NODE); 
    dof = 2;        sdof = n_node*dof;
    
    
    
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
    x_middle_up_up = abs(NODE(:, 1) - max(NODE(:, 1)) * 7/8) < 1e-6;  %% points of highest and lowest x,y %% points of highest and lowest x,y
    

    
    BC = zeros(n_node, 2);
    BC_N = zeros(n_node, 2);
    BC(x_lower,:) = 1;
    BCid = find(reshape(BC', [],1));
    freedofs = find(reshape(~BC', [],1));
    %% Force condition %%
    BC_N(x_upper & y_lower,1) = 1;
%     BC_N(x_middle & y_middle_low,1) = 1;
%     BC_N(x_middle_up & y_middle,2) = -1;

    BC_Nid = find(reshape(BC_N', [],1));
    p = 0.01; h = 25; l = 250;
    F = sparse(BC_Nid,1,1,sdof,1);
    
    F(BC_Nid(:)) = 10e3;                                                %p*(l/(n-1)) * t;
    
    %%%%%%%%%%%%%%%%%%%% Solve Stiffness %%%%%%%%%%%%%%%%%%%%
    E = 200e3; v = 0.25;
    D = E/(1-v^2) * [1 v 0; v 1 0; 0 0 (1-v)/2];
    [KE,edofMat] = stiffness_compl(NODE, ELEM, D, h); % K = global stiffness, KE = local stiffness
    iK = reshape(kron(edofMat,ones(8,1))',64*nele,1);
    jK = reshape(kron(edofMat,ones(1,8))',64*nele,1);
    sK = zeros(64*nele,1);
    for i = 1:nele
        sK(64*(i-1)+1:64*i) = reshape(KE(:,:,i)*x(i)'.^penal,64,1);
    end
    K = sparse(iK,jK,sK); %K = (K+K')/2;
    U = zeros(sdof,1); % displacement -> 2n-1 : u vector / 2n : v vector
    U(freedofs) = K(freedofs,freedofs)\F(freedofs);
    
    
    %%%%%%%%%%%%%%%%%%%% Objective Function & Sensitivity Analysis %%%%%%%%%%%%%%%%%%%%
    Ue = zeros(8,1);
    c_compl = 0;
    for i=1:nele
        for j = 1:4
            Ue(2*j-1,1) = U(2*ELEM(i,j) - 1);
            Ue(2*j,1) = U(2*ELEM(i,j));
        end
        ce = Ue' * KE(:,:,i) * Ue;
        c_compl = c_compl + (1e-9+x(i)^penal*(1-1e-9)) * ce;
        dc_compl(i) = -penal*(1e-9+x(i)^(penal-1)*(1-1e-9)) *ce;
    end
