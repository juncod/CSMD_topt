function [f0val,df0dx,fval,dfdx]=Stress_minimize(NODE,ELEM,x,Hs,H,pl,q,p,volfrac)
nele = length(ELEM);
x(:)=(H*x(:))./Hs;
[pnorm,pnorm_sen,MISES]=Stress_3D_Sensitivity_Comp(x,NODE,ELEM,pl,q,p);
% Plot Density
figure(1); colormap(gray);
patch('Faces',ELEM,'Vertices',NODE,'FaceVertexCData',ceil(x.*-2),'FaceColor','flat','LineStyle','none'); axis equal; axis tight; axis off;
dv = ones(nele,1)/(nele); 
sen(:) = H*(pnorm_sen(:)./Hs);
dv(:) = H*(dv(:)./Hs);%.*min(abs(pnorm_sen))./min(min(abs(dv)));
fval=[mean(x(:))-volfrac];
dfdx=[dv(:)'];
df0dx=sen';
f0val=pnorm;