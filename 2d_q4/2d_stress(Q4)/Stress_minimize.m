function [f0val,df0dx,fval,dfdx,MISES_MAX]=Stress_minimize(NODE,ELEM,x,Hs,H,pl,q,p,volfrac)
nele = length(ELEM);
x(:)=(H*x(:))./Hs;
[pnorm,pnorm_sen,MISES]=Stress_3D_Sensitivity_Comp(x,NODE,ELEM,pl,q,p);
% Plot Density
plot_x = max(0,x(:)-0.1);
figure(1); colormap(gray);
patch('Faces',ELEM,'Vertices',NODE,'FaceVertexCData',-x,'FaceColor','flat','LineStyle','none'); axis equal; axis tight; axis off;
%ceil(x.*-2) or -x
dv = ones(nele,1)/(nele); 
sen(:) = H*(pnorm_sen(:)./Hs);
dv(:) = H*(dv(:)./Hs);
% sen(:) = H*(x(:).*pnorm_sen(:))./Hs./max(1e-3,x(:)); % => sensitivity filter
fval=[mean(x(:))-volfrac];
dfdx=[dv(:)'];
df0dx=sen';
f0val=pnorm;
MISES_MAX = max(MISES);