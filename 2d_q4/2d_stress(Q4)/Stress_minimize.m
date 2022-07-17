function [f0val,df0dx,fval,dfdx,MISES_MAX]=Stress_minimize(NODE,ELEM,x,Hs,H,pl,q,p,volfrac)
nele = length(ELEM);
x(:)=(H*x(:))./Hs;
[pnorm,pnorm_sen,MISES]=Stress_3D_Sensitivity_Comp(x,NODE,ELEM,pl,q,p);

% Plot Density
plot_x = -ceil(max(0,x(:)-0.1));
figure(1); colormap(gray);
patch('Faces',ELEM,'Vertices',NODE,'FaceVertexCData',-x,'FaceColor','flat','LineStyle','none'); axis equal; axis tight; axis off; caxis([-1 0]);
figure(2); colormap(jet);
patch('Faces',ELEM,'Vertices',NODE,'FaceVertexCData',MISES,'FaceColor','flat','LineStyle','none'); axis equal; axis tight; axis off; drawnow; colorbar; caxis([min(MISES) max(MISES)]);
 
sen(:) = (H*pnorm_sen(:))./Hs;
sen(:) = H*(x(:).*sen(:))./Hs./max(1e-3,x(:));
dv = ones(nele,1)/(nele);
dv(:) = (H*dv(:))./Hs;
fval=[mean(x(:))-volfrac];
dfdx=[dv(:)'];
df0dx=sen';
f0val=pnorm;
MISES_MAX = max(MISES);