function J = jacobian_(NODE, ELEM,s,t)
    J = zeros(4,1);
    node = zeros(4,2);
    node(:,1) = NODE(ELEM(1,:),1);
    node(:,2) = NODE(ELEM(1,:),2);
    for i =1:4
        a = (1/4) * (-(1-t(i)) * node(1,1) + (1-t(i)) * node(2,1) + (1+t(i)) * node(3,1) - (1+t(i)) * node(4,1));  
        b = (1/4) * (-(1-t(i)) * node(1,2) + (1-t(i)) * node(2,2) + (1+t(i)) * node(3,2) - (1+t(i)) * node(4,2));
        c = (1/4) * (-(1-s(i)) * node(1,1) - (1+s(i)) * node(2,1) + (1+s(i)) * node(3,1) + (1-s(i)) * node(4,1));
        d = (1/4) * (-(1-s(i)) * node(1,2) - (1+s(i)) * node(2,2) + (1+s(i)) * node(3,2) + (1-s(i)) * node(4,2));
        jacobian = [a,b;c,d;];
        J(i,1) = det(jacobian);
    end
end