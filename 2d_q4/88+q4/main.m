clear; clc; close all;
addpath('data');
[NODE,ELEM] = inp_('main.inp');
volfrac = 0.5; penal = 3; rmin = 0.3;
saveFileName = '2';
plotCut = 0.1;
x = topology(NODE,ELEM,volfrac,penal,rmin,saveFileName,plotCut);
%% Function
function x = topology(NODE,ELEM,volfrac,penal,rmin,saveFileName,plotCut)
    nele = length(ELEM);
    Ue = zeros(8,1);
    x(1:nele) = volfrac;
    iter = 0;
    maxiter = 120;
    change = 120;
    [Hs,H]=prepare_filter(rmin,NODE,ELEM);
    f1 = figure; colormap(gray); axis equal; axis tight; axis off; caxis([-1 0]);
    f3 = figure; colormap(gray); axis equal; axis tight; axis off; caxis([-1 0]);


    ELEM_dis = zeros(nele,2);
    for i = 1:nele
        ELEM_dis(i,:) = (NODE(ELEM(i,1),:)+NODE(ELEM(i,2),:)+NODE(ELEM(i,3),:)+NODE(ELEM(i,4),:))/4;
    end
    eX_middle_upper = ELEM_dis(:,1)>5;
    eY_middle_upper = ELEM_dis(:,2)>5; 
    eZero = zeros(nele,1);
    eZero(eX_middle_upper&eY_middle_upper) = 1;
    eZero_id = find(reshape(eZero', [],1));

    % Start Iteration
    while change > 0.001 && iter < maxiter
        iter = iter + 1;
        if iter <= 15, gf = 1; else gf = min(2,1.01*gf); end
        x(eZero_id) = 0.001;
        xold = x;
    % FEA analysis
        [U,KE] = FE_(NODE,ELEM,x,penal);
    % Objective Function & Sensitivity Analysis
        c = 0;
        for i=1:nele
            for j = 1:4
                Ue(2*j-1,1) = U(2*ELEM(i,j) - 1);
                Ue(2*j,1) = U(2*ELEM(i,j));
            end
            c = c + (x(i)^penal) * Ue' * KE(:,:,i) * Ue;
            dc(i) = -penal*x(i)^(penal-1) *Ue' * KE(:,:,i) * Ue;
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
        plot_x = -ceil(max(0,x(:)-plotCut));
        figure(f1); 
        patch('Faces',ELEM,'Vertices',NODE,'FaceVertexCData',-x','FaceColor','flat','LineStyle','none'); 
        figure(f3); 
        patch('Faces',ELEM,'Vertices',NODE,'FaceVertexCData',plot_x,'FaceColor','flat','LineStyle','none'); 
        pause(1e-6);       
    end
    saveFunction(saveFileName,x,f1,f3,plotCut)
end


function [U,KE] = FE_(NODE,ELEM,x,penal)
    n_node = length(NODE); 
    nele = length(ELEM);
    dof = 2;        sdof = n_node*dof;
    %% Boundary Condition %%
    x_lower = abs(NODE(:,1)) < 1e-6;
    x_upper = abs(NODE(:,1)-10) < 1e-6;
    y_lower = abs(NODE(:,2)) < 1e-6;
    y_upper = abs(NODE(:,2)-10) < 1e-6;
    x_middle = abs(NODE(:,1)-5) < 1e-6;
    y_middle = abs(NODE(:,2)-5) < 1e-6;  %% points of highest and lowest x,y %% points of highest and lowest x,y

    BC = zeros(n_node, 2);
    BC(y_upper,:) = 1;
    BCid = find(reshape(BC', [],1));
    freedofs = find(reshape(~BC', [],1));
    %% Force condition %%
    BC_N = zeros(n_node, 2);
    BC_N(x_upper & y_middle,2) = 1;
    BC_Nid = find(reshape(BC_N', [],1));
    p = 0.01; h = 25; l = 250;
    F = sparse(BC_Nid,1,1,sdof,1);

    F(BC_Nid(:)) = -10e3;                                                %p*(l/(n-1)) * t;

    %% Solve Stiffness %%
    E = 200e3; v = 0.25;
    D = E/(1-v^2) * [1 v 0; v 1 0; 0 0 (1-v)/2];
    [KE,edofMat] = stiffness_(NODE, ELEM, D, h); % K = global stiffness, KE = local stiffness
    iK = reshape(kron(edofMat,ones(8,1))',64*nele,1);
    jK = reshape(kron(edofMat,ones(1,8))',64*nele,1);
    sK = zeros(64*nele,1);
    for i = 1:nele
        sK(64*(i-1)+1:64*i) = reshape(KE(:,:,i)*x(i)'.^penal,64,1);
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
        it = intersect(ex,ey);
        for j = 1:length(it)
            k=k+1;
            A = [m_elem(i,1), m_elem(i,2)];
            B = [m_elem(it(j),1), m_elem(it(j),2)];
            iH(k) = i;
            jH(k) = it(j);
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
