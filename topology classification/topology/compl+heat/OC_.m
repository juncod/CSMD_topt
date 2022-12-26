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
