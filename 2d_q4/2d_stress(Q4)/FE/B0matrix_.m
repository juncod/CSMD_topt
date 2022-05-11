function B0 = B0matrix_(NODE,ELEM,J0)
    %%% B0 matrix of node per element
    s=0;t=0;
    node = zeros(4,2);
    B0 = zeros(3,8);
    N=zeros(4,2);
    for i = 1:4
        node(i,1) = NODE(ELEM(1,i),1);
        node(i,2) = NODE(ELEM(1,i),2);
    end  
    N(1,1) = (t-1)/4;  N(1,2) = (s-1)/4;
    N(2,1) = (1-t)/4;  N(2,2) = -(1+s)/4;
    N(3,1) = (t+1)/4; N(3,2) = (s+1)/4;
    N(4,1) = -(t+1)/4; N(4,2) = (1-s)/4;
    a = (1/4) * (node(1,2)*(s-1)+node(2,2)*(-1-s)+node(3,2)*(1+s)+node(4,2)*(1-s));
    b = (1/4) * (node(1,2)*(t-1)+node(2,2)*(1-t)+node(3,2)*(1+t)+node(4,2)*(-1-t));
    c = (1/4) * (node(1,1)*(t-1)+node(2,1)*(1-t)+node(3,1)*(1+t)+node(4,1)*(-1-t));
    d = (1/4) * (node(1,1)*(s-1)+node(2,1)*(-1-s)+node(3,1)*(1+s)+node(4,1)*(1-s));
    for i=1:4      
        B0(1,2*i-1) = a*N(i,1) - b*N(i,2);
        B0(2,2*i) = c*N(i,2) - d*N(i,1);
        B0(3,2*i-1) = B0(2,2*i);
        B0(3,2*i) = B0(1,2*i-1);
    end
    B0 = B0 / J0;


end