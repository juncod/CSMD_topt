function [sen_H]=Sens_filter(NODE,ELEM,x)

    nele = length(ELEM);
    for i = 1:nele
        i_elem=ELEM(i);
        where_elem = ELEM==i_elem;
        neighbor_elem = find(sum(where_elem,2));
        sen_H(i)=sum(x(neighbor_elem))/length(neighbor_elem);
    end