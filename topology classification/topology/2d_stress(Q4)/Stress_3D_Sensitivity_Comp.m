function [pnorm,pnorm_sen,MISES,c,dc]=Stress_3D_Sensitivity_Comp(x,NODE,ELEM,pl,q,p)
n_node = length(NODE); nele = length(ELEM);
dof = 2;        sdof = n_node*dof;
%% Boundary Condition %%
x_lower = abs(NODE(:,1)) < 1e-6;
x_upper = abs(NODE(:,1)-10) < 1e-6;
y_lower = abs(NODE(:,2)) < 1e-6;
y_upper = abs(NODE(:,2)-10) < 1e-6;
y_middle = abs(NODE(:,2)-5) < 1e-6; %% points of highest and lowest x,y

BC = zeros(n_node, 2);
BC(y_upper,:) = 1;
BCid = find(reshape(BC', [],1));
freedofs = find(reshape(~BC', [],1));
%% Force condition %%
BC_N = zeros(n_node, 2);
BC_N(x_upper & y_middle,2) = 1;
BC_Nid = find(reshape(BC_N', [],1));
h = 10;
F = sparse(BC_Nid,1,1,sdof,1);
F(BC_Nid(:)) = -10e3;
U = zeros(sdof,1);
%% Solve Stiffness %%
E0 = 200e3; Emin =1e-9; v = 0.25;
D = 1/(1-v^2) * [1 v 0; v 1 0; 0 0 (1-v)/2];
[KE,edofMat,B0] = stiffness_(NODE, ELEM, D, h); % K = global stiffness, KE = local stiffness
iK = reshape(kron(edofMat,ones(8,1))',64*nele,1);
jK = reshape(kron(edofMat,ones(1,8))',64*nele,1);
sK = zeros(64*nele,1);
for i = 1:nele
    sK(64*(i-1)+1:64*i) = reshape(KE(:,:,i)*(Emin+x(i)'.^pl*(E0-Emin)),64,1);
end
K = sparse(iK,jK,sK); %K = (K+K')/2;
U(freedofs) = K(freedofs,freedofs)\F(freedofs);
% U", K", KE(3차원), edofMat", D" , B(B0), freedofs"


Ue = zeros(8,1);
c = 0;
for i=1:nele
    for j = 1:4
        Ue(2*j-1,1) = U(2*ELEM(i,j) - 1);
        Ue(2*j,1) = U(2*ELEM(i,j));
    end
    c = c + (Emin+x(i)'.^pl*(E0-Emin)) * Ue' * KE(:,:,i) * Ue;
    dc(i) = -pl*x(i)^(pl-1)*(E0-Emin) *Ue' * KE(:,:,i) * Ue;
end




MISES=zeros(nele,1); % von Mises stress vector
S = zeros(nele,3);
mises_mat = [1 -1/2 0; -1/2 1 0; 0 0 3];
for i=1:nele
    S(i,:)=x(i)^q*((Emin+x(i)'.^pl*(E0-Emin))*D*B0(:,:,i)*U(edofMat(i,:)));
    MISES(i)=sqrt(S(i,:)*mises_mat*S(i,:)');
end
index_matrix=edofMat';
pnorm=(sum(MISES.^p))^(1/p);




gama=zeros(sdof,1);
for i=1:nele
    index=index_matrix(:,i);
    gama(index)=gama(index)+x(i)^q*B0(:,:,i)'*D'*mises_mat*S(i,:)'*MISES(i).^(p-2);
end
lamda=zeros(sdof,1);
lamda(freedofs,:)=K(freedofs,freedofs)\gama(freedofs,:);


pnorm_sen = zeros(nele,1);
for i = 1:nele
    index=index_matrix(:,i);
    pnorm_sen(i) = (pnorm^(1-p)*(MISES(i)^p*q/x(i)-pl*x(i)^(pl-1)*lamda(index)'*KE(:,:,i)*U(index)));
end



% =========================================================================
% The present code is part of the journal paper by Deng et al. 2021 and is
% derived from the code which was part of the paper by Liu et al. 2014.
% -------------------------------------------------------------------------
% Please send your suggestions and comments to: albertto@pitt.edu
% -------------------------------------------------------------------------
% The code is intended for educations purposes, and the details and
% extensions can be found in the paper:
% Deng, H., Vulimiri, P.S. & To, A.C. An efficient 146-line 3D sensitivity
% analysis code of stress-based topology optimization written in MATLAB.
% Optim Eng (2021). https://doi.org/10.1007/s11081-021-09675-3
% -------------------------------------------------------------------------
% Details of the finite element formulation from Liu et al. can be found in
% the paper:
% Liu, K., Tovar, A. An efficient 3D topology optimization code written in
% Matlab. Struct Multidisc Optim 50, 1175–1196 (2014).
% https://doi.org/10.1007/s00158-014-1107-x
% -------------------------------------------------------------------------
% The code can be downloaded from the website:
% https://github.com/PittAMRL/StressTopOpt
% -------------------------------------------------------------------------
% The code from Liu et al. from which this is derived from, can be
% downloaded from: http://www.top3dapp.com/
% -------------------------------------------------------------------------
% Disclaimer:
% The code may be distributed and used for educational purposes.
% The authors do not guarantee that the code is free from errors.
