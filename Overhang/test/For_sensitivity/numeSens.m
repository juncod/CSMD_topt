function [ ndc ] = numeSens( x ,penal)
    startX=x;
    [nely,nelx,nelz] = size(x);
    
    % USER-DEFINED MATERIAL PROPERTIES
    E0 = 1;           % Young's modulus of solid material
    Emin = 1e-9;      % Young's modulus of void-like material
    nu = 0.3;         % Poisson's ratio
    % USER-DEFINED LOAD DOFs
    [il,jl,kl] = meshgrid(nelx, 0, 0:nelz);                 % Coordinates
    loadnid = kl*(nelx+1)*(nely+1)+il*(nely+1)+(nely+1-jl); % Node IDs
    loaddof = 3*loadnid(:) - 1;                             % y방향 DOF
    % USER-DEFINED SUPPORT FIXED DOFs
    [iif,jf,kf] = meshgrid(0,0:nely,0:nelz);                  % Coordinates
    fixednid = kf*(nelx+1)*(nely+1)+iif*(nely+1)+(nely+1-jf); % Node IDs
    fixeddof = [3*fixednid(:); 3*fixednid(:)-1; 3*fixednid(:)-2]; % xyz방향 DOF
    % PREPARE FINITE ELEMENT ANALYSIS
    nele = nelx*nely*nelz;
    ndof = 3*(nelx+1)*(nely+1)*(nelz+1);
    F = sparse(loaddof,1,-1,ndof,1); % ==> y방향 아래로 힘
    U = zeros(ndof,1);
    freedofs = setdiff(1:ndof,fixeddof);
    KE = lk_H8(nu);
    nodegrd = reshape(1:(nely+1)*(nelx+1),nely+1,nelx+1); % ==> 한 층의 node ids
    nodeids = reshape(nodegrd(1:end-1,1:end-1),nely*nelx,1); % ==> 한 층의 element당 첫번째 node ids 
    nodeidz = 0:(nely+1)*(nelx+1):(nelz-1)*(nely+1)*(nelx+1); % ==> 층끼리 node id간의 차이
    nodeids = repmat(nodeids,size(nodeidz))+repmat(nodeidz,size(nodeids)); % ==> 모든 element당 첫번째 node ids
    edofVec = 3*nodeids(:)+1; % ==> 모든 element당 두번째 node dof_x
    edofMat = repmat(edofVec,1,24)+ ... % ==> 모든 element당 모든 node의 dofs
        repmat([0 1 2 3*nely + [3 4 5 0 1 2] -3 -2 -1 ...
        3*(nely+1)*(nelx+1)+[0 1 2 3*nely + [3 4 5 0 1 2] -3 -2 -1]],nele,1);
    iK = reshape(kron(edofMat,ones(24,1))',24*24*nele,1);
    jK = reshape(kron(edofMat,ones(1,24))',24*24*nele,1);

    % FE-ANALYSIS
    sK = reshape(KE(:)*(Emin+x(:)'.^penal*(E0-Emin)),24*24*nele,1); % ==> (15)
    K = sparse(iK,jK,sK); K = (K+K')/2;
    U(freedofs,:) = K(freedofs,freedofs)\F(freedofs,:); % ==> U = inv(K)*F
    % OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
    ce = reshape(sum((U(edofMat)*KE).*U(edofMat),2),[nely,nelx,nelz]); % ==> UKU
    ci = sum(sum(sum((Emin+x.^penal*(E0-Emin)).*ce))); % ==> c=E(x)*ce


    ndc=zeros([nely,nelx,nelz]);
    for i=1:nely
    for j=1:nelx
    for k=1:nelz
        x=startX;
        x(i,j,k) = x(i,j,k) + 0.0001;
        % USER-DEFINED MATERIAL PROPERTIES
        E0 = 1;           % Young's modulus of solid material
        Emin = 1e-9;      % Young's modulus of void-like material
        nu = 0.3;         % Poisson's ratio
        % USER-DEFINED LOAD DOFs
        [il,jl,kl] = meshgrid(nelx, 0, 0:nelz);                 % Coordinates
        loadnid = kl*(nelx+1)*(nely+1)+il*(nely+1)+(nely+1-jl); % Node IDs
        loaddof = 3*loadnid(:) - 1;                             % y방향 DOF
        % USER-DEFINED SUPPORT FIXED DOFs
        [iif,jf,kf] = meshgrid(0,0:nely,0:nelz);                  % Coordinates
        fixednid = kf*(nelx+1)*(nely+1)+iif*(nely+1)+(nely+1-jf); % Node IDs
        fixeddof = [3*fixednid(:); 3*fixednid(:)-1; 3*fixednid(:)-2]; % xyz방향 DOF
        % PREPARE FINITE ELEMENT ANALYSIS
        nele = nelx*nely*nelz;
        ndof = 3*(nelx+1)*(nely+1)*(nelz+1);
        F = sparse(loaddof,1,-1,ndof,1); % ==> y방향 아래로 힘
        U = zeros(ndof,1);
        freedofs = setdiff(1:ndof,fixeddof);
        KE = lk_H8(nu);
        nodegrd = reshape(1:(nely+1)*(nelx+1),nely+1,nelx+1); % ==> 한 층의 node ids
        nodeids = reshape(nodegrd(1:end-1,1:end-1),nely*nelx,1); % ==> 한 층의 element당 첫번째 node ids 
        nodeidz = 0:(nely+1)*(nelx+1):(nelz-1)*(nely+1)*(nelx+1); % ==> 층끼리 node id간의 차이
        nodeids = repmat(nodeids,size(nodeidz))+repmat(nodeidz,size(nodeids)); % ==> 모든 element당 첫번째 node ids
        edofVec = 3*nodeids(:)+1; % ==> 모든 element당 두번째 node dof_x
        edofMat = repmat(edofVec,1,24)+ ... % ==> 모든 element당 모든 node의 dofs
            repmat([0 1 2 3*nely + [3 4 5 0 1 2] -3 -2 -1 ...
            3*(nely+1)*(nelx+1)+[0 1 2 3*nely + [3 4 5 0 1 2] -3 -2 -1]],nele,1);

        % FE-ANALYSIS
        sK = reshape(KE(:)*(Emin+x(:)'.^penal*(E0-Emin)),24*24*nele,1); % ==> (15)
        K = sparse(iK,jK,sK); K = (K+K')/2;
        U(freedofs,:) = K(freedofs,freedofs)\F(freedofs,:); % ==> U = inv(K)*F
        % OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
        ce = reshape(sum((U(edofMat)*KE).*U(edofMat),2),[nely,nelx,nelz]); % ==> UKU
        c = sum(sum(sum((Emin+x.^penal*(E0-Emin)).*ce))); % ==> c=E(x)*ce
        ndc(i,j,k) = (c-ci)/(0.0001);
    end
    end
    end


end


% === GENERATE ELEMENT STIFFNESS MATRIX ===
function [KE] = lk_H8(nu)
    A = [32 6 -8 6 -6 4 3 -6 -10 3 -3 -3 -4 -8;
        -48 0 0 -24 24 0 0 0 12 -12 0 12 12 12];
    k = 1/144*A'*[1; nu];
    
    K1 = [k(1) k(2) k(2) k(3) k(5) k(5);
        k(2) k(1) k(2) k(4) k(6) k(7);
        k(2) k(2) k(1) k(4) k(7) k(6);
        k(3) k(4) k(4) k(1) k(8) k(8);
        k(5) k(6) k(7) k(8) k(1) k(2);
        k(5) k(7) k(6) k(8) k(2) k(1)];
    K2 = [k(9)  k(8)  k(12) k(6)  k(4)  k(7);
        k(8)  k(9)  k(12) k(5)  k(3)  k(5);
        k(10) k(10) k(13) k(7)  k(4)  k(6);
        k(6)  k(5)  k(11) k(9)  k(2)  k(10);
        k(4)  k(3)  k(5)  k(2)  k(9)  k(12)
        k(11) k(4)  k(6)  k(12) k(10) k(13)];
    K3 = [k(6)  k(7)  k(4)  k(9)  k(12) k(8);
        k(7)  k(6)  k(4)  k(10) k(13) k(10);
        k(5)  k(5)  k(3)  k(8)  k(12) k(9);
        k(9)  k(10) k(2)  k(6)  k(11) k(5);
        k(12) k(13) k(10) k(11) k(6)  k(4);
        k(2)  k(12) k(9)  k(4)  k(5)  k(3)];
    K4 = [k(14) k(11) k(11) k(13) k(10) k(10);
        k(11) k(14) k(11) k(12) k(9)  k(8);
        k(11) k(11) k(14) k(12) k(8)  k(9);
        k(13) k(12) k(12) k(14) k(7)  k(7);
        k(10) k(9)  k(8)  k(7)  k(14) k(11);
        k(10) k(8)  k(9)  k(7)  k(11) k(14)];
    K5 = [k(1) k(2)  k(8)  k(3) k(5)  k(4);
        k(2) k(1)  k(8)  k(4) k(6)  k(11);
        k(8) k(8)  k(1)  k(5) k(11) k(6);
        k(3) k(4)  k(5)  k(1) k(8)  k(2);
        k(5) k(6)  k(11) k(8) k(1)  k(8);
        k(4) k(11) k(6)  k(2) k(8)  k(1)];
    K6 = [k(14) k(11) k(7)  k(13) k(10) k(12);
        k(11) k(14) k(7)  k(12) k(9)  k(2);
        k(7)  k(7)  k(14) k(10) k(2)  k(9);
        k(13) k(12) k(10) k(14) k(7)  k(11);
        k(10) k(9)  k(2)  k(7)  k(14) k(7);
        k(12) k(2)  k(9)  k(11) k(7)  k(14)];
    KE = 1/((nu+1)*(1-2*nu))*...
        [ K1  K2  K3  K4;
        K2'  K5  K6  K3';
        K3' K6  K5' K2';
        K4  K3  K2  K1'];
    end