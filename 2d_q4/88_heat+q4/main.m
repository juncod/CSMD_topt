clear; clc;
tic
[NODE,ELEM] = inp_('Job-2.inp');
volfrac = 0.5; penal = 3; rmin = 0.3;
x = topology(NODE,ELEM,volfrac,penal,rmin);
toc
%% Function
function x = topology(NODE,ELEM,volfrac,penal,rmin)
    k0 = 1; kmin = 1e-3;
    nele = length(ELEM);
    Ue = zeros(4,1);
    x(1:nele) = volfrac;
    iter = 0;
    maxiter = 150;
    change = 1;
    [Hs,H]=prepare_filter(rmin,NODE,ELEM);
    % Start Iteration
    while change > 0.001 && iter < maxiter
        iter = iter + 1;
        if iter <= 15, gf = 1; else gf = min(1.5,1.01*gf); end
        xold = x;
    % FEA analysis
        [U,KE] = FE_(NODE,ELEM,x,penal,k0,kmin);
    % Objective Function & Sensitivity Analysis
        c = 0;
        for i=1:nele
            for j = 1:4
                Ue(j,1) = U(ELEM(i,j));
            end
            c = c + (kmin+x(i)^penal*(k0-kmin)) * Ue' * KE(:,:,i) * Ue;
            dc(i) = -penal*x(i)^(penal-1)*(k0-kmin) *Ue' * KE(:,:,i) * Ue;
        end
        % Filtering of Sensitivities
        dcn(:) = H*(x(:).*dc(:))./Hs./max(1e-3,x(:));
        % Design Update by the Optimality Criteria Method
        [x] = OC_(ELEM,x,volfrac,dcn,gf);
        % Print Results
        change = max(max(abs(x-xold)));
        disp([' It.: ' sprintf('%4i',iter) ' Obj.: ' sprintf('%10.4f',c) ...
       ' Vol.: ' sprintf('%6.3f',sum(sum(x))/(nele)) ...
        ' ch.: ' sprintf('%6.3f',change)])
        % Plot Density
        patch('Faces',ELEM,'Vertices',NODE,'FaceVertexCData',-x','FaceColor','flat','LineStyle','none'); axis equal; axis tight; axis off;
        colormap(gray)
        pause(1e-6);       
    end
end


function [U,KE] = FE_(NODE,ELEM,x,penal,k0,kmin)
    n_node = length(NODE); 
    nele = length(ELEM);
    dof = 1;        sdof = n_node*dof;
    %% Boundary Condition %%
    x_lower = abs(NODE(:,1)) < 1e-6;
    x_upper = abs(NODE(:,1)-max(NODE(:,1))) < 1e-6;
    y_lower = abs(NODE(:,2)) < 1e-6;
    y_upper = abs(NODE(:,2)-max(NODE(:,2))) < 1e-6;
    x_middle = abs(NODE(:,1)-max(NODE(:,1)/2)) < 0.5;
    y_middle = abs(NODE(:,2)-max(NODE(:,2)/2)) < 0.5;
    BC_nodes = y_upper & x_middle;
    % BC_nodes = x_upper & y_lower;

    BC = zeros(n_node, 1);
    BC(BC_nodes,1) = 1;
    BCid = find(reshape(BC', [],1));
    freedofs = find(reshape(~BC', [],1));

    BC_N = zeros(n_node, 1);
    BC_N(x_lower,1) = 1;
    BC_Nid = find(reshape(BC_N', [],1));
    %% Force condition %%
    F = zeros(sdof,1);% => q, 열이 빠지는 곳
    % F(BC_Nid(:)) = -0.01;
    F(:) = -0.01; 
    %% Solve Stiffness %%
    v = 0.25;    h = 25;
    D = 1/(1-v^2) * [1 v 0; v 1 0; 0 0 (1-v)/2];
    [KE,edofMat] = stiffness_(NODE, ELEM, D, h); % K = global stiffness, KE = local stiffness
    iK = reshape(kron(edofMat,ones(4,1))',16*nele,1);
    jK = reshape(kron(edofMat,ones(1,4))',16*nele,1);
    sK = zeros(16*nele,1);
    for i = 1:nele
        sK(16*(i-1)+1:16*i) = reshape(KE(:,:,i)*(kmin+x(i)^penal*(k0-kmin)),16,1);
    end
    K = sparse(iK,jK,sK); %K = (K+K')/2;
    U = zeros(sdof,1); % displacement -> 2n-1 : u vector / 2n : v vector
    U(freedofs) = K(freedofs,freedofs)\F(freedofs);
end

%% mesh independcy filtering
function [Hs,H]=prepare_filter(rmin,NODE,ELEM)
    nele=length(ELEM);
    iH = zeros(nele*nele,1); % 그냥 최대치 설정
    jH = zeros(size(iH));
    sH = zeros(size(iH));
    m_elem = zeros(nele,2);
    for i = 1:nele
        m_elem(i,1) = m_elem(i,1) + sum(NODE(ELEM(i,:),1))/4;
        m_elem(i,2) = m_elem(i,2) + sum(NODE(ELEM(i,:),2))/4;
    end % => m_elem : element의 x좌표 y좌표
    k = 0;
    for i = 1:nele
        ex = find(abs(m_elem(:,1) - m_elem(i,1)) < rmin);
        ey = find(abs(m_elem(:,2) - m_elem(i,2)) < rmin);
        E = intersect(ex,ey);
        for j = 1:length(E)
            k=k+1;
            A = [m_elem(i,1), m_elem(i,2)];
            B = [m_elem(E(j),1), m_elem(E(j),2)];
            iH(k) = i;
            jH(k) = E(j);
            sH(k) = max(0,rmin - norm(A-B));
        end
    end
    H = sparse(iH(1:k),jH(1:k),sH(1:k));
    Hs = sum(H,2);
end




%% x new
function [xnew] = OC_(ELEM,x,volfrac,dcn,gf)
    l1 = 0; l2 = 1e5; move = 0.05;
    nele = length(ELEM);
    dv = ones(1,nele)/nele;
    while(l2-l1 > 1e-6)
        lmid = 0.5*(l1+l2);
        xnew = max(0.001,max(x-move,min(1,min(x+move,(x.*sqrt(-dcn./dv./lmid)).^gf))));
        if sum(sum(xnew)) - volfrac*nele > 0
            l1 = lmid;
        else
            l2 = lmid;
        end
    end
end
