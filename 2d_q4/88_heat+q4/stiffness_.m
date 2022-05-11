function [KE,edofMat] = stiffness_(NODE, ELEM, D, h)
k0   = 1;
nele = size(ELEM,1);
nnel = 4;                               % number of nodes per element
ndof = 1;                               % number of dofs per node
n_node = size(NODE,1);                  % total number of nodes
sdof = n_node*ndof;                      % total number of dofs
edofMat = zeros(nele,4);
root3 = 1/sqrt(3);
s = (root3)*[-1,1,1,-1];
t = (root3)*[-1,-1,1,1];
KE = zeros(4,4,nele);
for ii = 1:nele
    J = jacobian_(NODE,ELEM(ii,:),s,t);
    Ke = zeros(4,4);
    for i = 1:4
        B = Bmatrix_(NODE,ELEM(ii,:),J(i,1),s(i),t(i));
        Ke(:,:) = Ke(:,:) + B' * [k0 0; 0 k0] * B * J(i,1) * h;        
    end
    KE(:,:,ii) = Ke(:,:);

    nodeId = zeros(1,4);
    dofId = zeros(1,4);
    nodeId(:) = ELEM(ii,:);
    dofId(:) = [nodeId(1) nodeId(2) nodeId(3) nodeId(4)];
    edofMat(ii,:) = dofId;

    %%
end

