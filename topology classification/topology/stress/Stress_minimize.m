function [f0val,df0dx,fval,dfdx,MISES,sen,x]=Stress_minimize(NODE,ELEM,x,Hs,H,pl,q,p,volfrac)
nele = length(ELEM);
x(:)=(H*x(:))./Hs;
[pnorm,pnorm_sen,MISES]=Stress_3D_Sensitivity_Comp(x,NODE,ELEM,pl,q,p);

% x(ex_id(:))=1;
% pnorm_sen(ex_id(:))=0;
sen(:) = (H*pnorm_sen(:))./Hs;
% sen(ex_id(:))=0;
% sen(:) = H*(x(:).*sen(:))./Hs./max(1e-3,x(:));
dv = ones(nele,1)/(nele);
dv(:) = (H*dv(:))./Hs;


fval=[mean(x(:))-volfrac];
dfdx=[dv(:)'];
df0dx=sen';
f0val=pnorm;