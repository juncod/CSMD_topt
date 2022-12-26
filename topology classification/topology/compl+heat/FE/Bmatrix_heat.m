function B = Bmatrix_heat(NODE,ELEM,J,s,t)
%%% B matrix of node per element
node = zeros(4,2);
B = zeros(2,4);
N=zeros(4,2);
N(1,1) = (t-1)/4;  N(1,2) = (s-1)/4;
N(2,1) = (1-t)/4;  N(2,2) = -(1+s)/4;
N(3,1) = (t+1)/4; N(3,2) = (s+1)/4;
N(4,1) = -(t+1)/4; N(4,2) = (1-s)/4;
node(:,1) = NODE(ELEM(1,:),1);
node(:,2) = NODE(ELEM(1,:),2); 
a = (1/4) * (node(1,2)*(s-1)+node(2,2)*(-1-s)+node(3,2)*(1+s)+node(4,2)*(1-s));
b = (1/4) * (node(1,2)*(t-1)+node(2,2)*(1-t)+node(3,2)*(1+t)+node(4,2)*(-1-t));
c = (1/4) * (node(1,1)*(t-1)+node(2,1)*(1-t)+node(3,1)*(1+t)+node(4,1)*(-1-t));
d = (1/4) * (node(1,1)*(s-1)+node(2,1)*(-1-s)+node(3,1)*(1+s)+node(4,1)*(1-s));

for i=1:4      
    B(1,i) = a*N(i,1) - b*N(i,2);
    B(2,i) = c*N(i,2) - d*N(i,1);
end
B = B / J;
