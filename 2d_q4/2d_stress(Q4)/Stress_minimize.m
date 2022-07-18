function [f0val,df0dx,fval,dfdx,MISES,pnorm]=Stress_minimize(NODE,ELEM,x,Hs,H,pl,q,p,volfrac)
nele = length(ELEM);
% x(:)=(H*x(:))./Hs;
[pnorm,pnorm_sen,MISES,c,dc]=Stress_3D_Sensitivity_Comp(x,NODE,ELEM,pl,q,p);


% sen(:) = (H*pnorm_sen(:))./Hs;
% dcn(:) = (H*dc(:))./Hs;
dcn(:) = H*(x(:).*dc(:))./Hs./max(1e-3,x(:));
dv = ones(nele,1)/(nele);
% dv(:) = (H*dv(:))./Hs;
maxStr = 5000;
fval=[mean(x(:))-volfrac; pnorm-maxStr];
dfdx=[dv(:)';pnorm_sen(:)'];
f0val=c;
df0dx=dcn';
