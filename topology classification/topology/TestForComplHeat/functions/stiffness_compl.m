function [KE,edofMat] = stiffness_compl(NODE, ELEM, D, h)
nele = size(ELEM,1);
nnel = 4;                               % number of nodes per element
ndof = 2;                               % number of dofs per node
n_node = size(NODE,1);            % total number of nodes
sdof = n_node*ndof;                      % total number of dofs
edofMat = zeros(nele,8);
root3 = 1/sqrt(3);
s = (root3)*[-1,1,1,-1];
t = (root3)*[-1,-1,1,1];
KE = zeros(8,8,nele);
for ii = 1:nele
    J = jacobian_compl(NODE,ELEM(ii,:),s,t);
    Ke = zeros(8,8);
    for i = 1:4
        B = Bmatrix_compl(NODE,ELEM(ii,:),J(i,1),s(i),t(i));
        Ke(:,:) = Ke(:,:) + B' * D * B * J(i,1) * h;        
    end
    KE(:,:,ii) = Ke(:,:);

    nodeId = zeros(nnel,1);
    dofId = zeros(nnel*ndof,1);
    nodeId(:) = ELEM(ii,:);
    dofId(:) = [nodeId(1)*ndof-1 nodeId(1)*ndof nodeId(2)*ndof-1 nodeId(2)*ndof ...
                    nodeId(3)*ndof-1 nodeId(3)*ndof nodeId(4)*ndof-1 nodeId(4)*ndof];
    edofMat(ii,:) = dofId;
end

