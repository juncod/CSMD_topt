function [KE,edofMat,B0] = stiffness_(NODE, ELEM, D, h)
nele = size(ELEM,1);
nnel = 4;                               % number of nodes per element
ndof = 2;                               % number of dofs per node
n_node = size(NODE,1);                  % total number of nodes
sdof = n_node*ndof;                     % total number of dofs

root3 = 1/sqrt(3);
s = (root3)*[-1,1,1,-1];
t = (root3)*[-1,-1,1,1];
KE = zeros(8,8,nele);
B0 = zeros(3,8,nele);
edofMat = zeros(nele,8);
for ii = 1:nele
    a = zeros(3,8);
    [J,J0] = jacobian_(NODE,ELEM(ii,:),s,t);
    Ke = zeros(8,8);
    for i = 1:4
        B = Bmatrix_(NODE,ELEM(ii,:),J(i,1),s(i),t(i));
        Ke(:,:) = Ke(:,:) + B' * D * B * J(i,1) * h;        
        a = a(:) + B(:);
    end
    a = a/4;
    B0(:,:,ii) = B0matrix_(NODE,ELEM(ii,:),J0);
    KE(:,:,ii) = Ke(:,:);

    nodeId = zeros(nnel,1);
    dofId = zeros(nnel*ndof,1);
    nodeId(:) = ELEM(ii,:);
    dofId(:) = [nodeId(1)*ndof-1 nodeId(1)*ndof nodeId(2)*ndof-1 nodeId(2)*ndof ...
                    nodeId(3)*ndof-1 nodeId(3)*ndof nodeId(4)*ndof-1 nodeId(4)*ndof];
    edofMat(ii,:) = dofId;
end