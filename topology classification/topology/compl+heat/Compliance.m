function [c_compl,dc_compl] = Compliance(x,NODE,ELEM,penal)

nele = length(ELEM);
n_node = length(NODE); 
dof = 2;        sdof = n_node*dof;



%%%%%%%%%%%%%%%%%%%% Boundary Condition %%%%%%%%%%%%%%%%%%%%
x_lower = abs(NODE(:,1)) < 1e-6;
x_upper = abs(NODE(:,1)-10) < 1e-6;
y_lower = abs(NODE(:,2)) < 1e-6;
y_upper = abs(NODE(:,2)-10) < 1e-6;
x_middle = abs(NODE(:,1)-5) < 1e-6;
y_middle = abs(NODE(:,2)-5) < 1e-6;  %% points of highest and lowest x,y %% points of highest and lowest x,y

x_point1 = 7.5; y_point1 = 10;
node_point1 = abs(NODE(:,1)-x_point1) < 1e-6 & abs(NODE(:,2)-y_point1) < 1e-6;
x_point2 = 7.5; y_point2 = 0;
node_point2 = abs(NODE(:,1)-x_point2) < 1e-6 & abs(NODE(:,2)-y_point2) < 1e-6;
x_point3 = 6; y_point3 = 5;
node_point3 = abs(NODE(:,1)-x_point3) < 1e-6 & abs(NODE(:,2)-y_point3) < 1e-6;
x_point4 = 6; y_point4 = 4;
node_point4 = abs(NODE(:,1)-x_point4) < 1e-6 & abs(NODE(:,2)-y_point4) < 1e-6;
x_point5 = 5; y_point5 = 5;
node_point5 = abs(NODE(:,1)-x_point5) < 1e-6 & abs(NODE(:,2)-y_point5) < 1e-6;

BC = zeros(n_node, 2);
BC(x_lower,:) = 1;
% BC(x_lower & y_lower ,:) = 1;
% BC(x_lower & y_upper ,:) = 1;
% BC(x_lower & y_middle ,:) = 1;
% BC(x_middle & y_lower ,:) = 1;
% BC(x_middle & y_upper ,:) = 1;
% BC(x_upper & y_lower ,:) = 1;
% BC(x_upper & y_middle ,:) = 1;
% BC(x_upper & y_upper ,:) = 1;
BCid = find(reshape(BC', [],1));
freedofs = find(reshape(~BC', [],1));
%% Force condition %%
BC_N(x_upper & y_middle,2) = -1;
BC_Nid = find(reshape(BC_N', [],1));
p = 0.01; h = 25; l = 250;
F = sparse(BC_Nid,1,1,sdof,1);

F(BC_Nid(:)) = -10e3;                                                %p*(l/(n-1)) * t;

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
