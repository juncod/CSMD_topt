function [Hs,H]=Prepare_filter(rmin,NODE,ELEM)
    nele=length(ELEM);
    iH = zeros(nele*nele,1); % for문 itreation 수
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