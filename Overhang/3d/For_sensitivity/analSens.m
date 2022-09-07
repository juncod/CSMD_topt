% AN 169 LINE 3D TOPOLOGY OPITMIZATION CODE BY LIU AND TOVAR (JUL 2013) 
% analSens(60,20,4,0.3,3,1.5,1)
function analSens(nelx,nely,nelz,volfrac,penal,rmin,qmax)
% USER-DEFINED LOOP PARAMETERS
maxloop = 200;    % Maximum number of iterations
tolx = 0.01;      % Terminarion criterion
displayflag = 0;  % Display structure flag
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
jK = reshape(kron(edofMat,ones(1,24))',24*24*nele,1); % ==> 모든 element의 모든 node의 dofs 24번 반복
% PREPARE FILTER
iH = ones(nele*(2*(ceil(rmin)-1)+1)^2,1);
jH = ones(size(iH));
sH = zeros(size(iH));
k = 0;
for k1 = 1:nelz
    for i1 = 1:nelx
        for j1 = 1:nely
            e1 = (k1-1)*nelx*nely + (i1-1)*nely+j1;
            for k2 = max(k1-(ceil(rmin)-1),1):min(k1+(ceil(rmin)-1),nelz)
                for i2 = max(i1-(ceil(rmin)-1),1):min(i1+(ceil(rmin)-1),nelx)
                    for j2 = max(j1-(ceil(rmin)-1),1):min(j1+(ceil(rmin)-1),nely)
                        e2 = (k2-1)*nelx*nely + (i2-1)*nely+j2;
                        k = k+1;
                        iH(k) = e1;
                        jH(k) = e2;
                        sH(k) = max(0,rmin-sqrt((i1-i2)^2+(j1-j2)^2+(k1-k2)^2));
                    end
                end
            end
        end
    end
end
H = sparse(iH,jH,sH);
Hs = sum(H,2);
% INITIALIZE ITERATION
x = repmat(volfrac,[nely,nelx,nelz]);
xPhys = x; baseplate='S';
[xPrint] = AM_filter(xPhys,baseplate);


% FE-ANALYSIS
sK = reshape(KE(:)*(Emin+xPrint(:)'.^penal*(E0-Emin)),24*24*nele,1); % ==> (15)
K = sparse(iK,jK,sK); K = (K+K')/2;
U(freedofs,:) = K(freedofs,freedofs)\F(freedofs,:); % ==> U = inv(K)*F
% OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
ce = reshape(sum((U(edofMat)*KE).*U(edofMat),2),[nely,nelx,nelz]); % ==> UKU
c = sum(sum(sum((Emin+xPrint.^penal*(E0-Emin)).*ce))); % ==> c=E(x)*ce
dc = -penal*(E0-Emin)*xPrint.^(penal-1).*ce; % ==> (27)
dv = ones(nely,nelx,nelz); % ==> (20)
% FILTERING AND MODIFICATION OF SENSITIVITIES
dc(:) = H*(dc(:)./Hs); % ==> (22)
dv(:) = H*(dv(:)./Hs); % ==> (19)
[xPrint, dc,dv] = AM_filter(xPhys,baseplate,dc,dv);
[ ndc ] = numeSens( x, penal);
anal = reshape(dc,1,nele);
nume = reshape(ndc,1,nele);
plot(1:nele,anal,'b-');hold on;plot(1:nele,nume,'r-');hold on;
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