clear; clc; close all;
addpath('FE'); addpath('MMA'); addpath('data');
% pl:penal, q:stress relaxation, p:p-norm
rmin = 0.5;   pl = 3;   q = 0.5;   p = 20;   volfrac = 0.5;

saveFileName = '4';
plotCut = 0.15;

[NODE,ELEM] = inp_('Job-L.inp');
[Hs,H]=Prepare_filter(rmin,NODE,ELEM);
nele = length(ELEM);
x=volfrac*ones(nele,1);
outeriter = 0;
f1 = figure; colormap(gray); axis equal; axis tight; axis off; caxis([-1 0]);
f2 = figure; colormap(jet); axis equal; axis tight; axis off; colorbar;
f3 = figure; colormap(gray); axis equal; axis tight; axis off; caxis([-1 0]);
%%%%%%%%%%%%%%%%%% M M A Zone %%%%%%%%%%%%%%%%%%
m =2;
epsimin = 0.0000001;
n=length(x(:));
xold1   = x;
xold2   = x;
xlb = 1e-3*ones(n,1);
xub = 1*ones(n,1);
xmin    = xlb;
xmax    = xub;
low     = xlb;
upp     = xub;
c       = [1e6 1e6]';
d       = [0 0]';
a0      = 1;
a       = [0 0]';
raa0    = 0.0001;
raa     = 0.0001;
raa0eps = 0.0000001;
raaeps  = 0.0000001;
maxoutit  = 150;
kkttol  = 0;
nele = length(ELEM);
x_his=zeros(nele,maxoutit);
kktnorm = kkttol+1;
outit = 0;
while  outit < maxoutit
    outit   = outit+1;
    if outit <= 15, gf = 0.2; else gf = min(0.5,1.01*gf); end
    outeriter = outeriter+1;
    [f0val,df0dx,fval,dfdx,MISES,pnorm]=Stress_minimize(NODE,ELEM,x,Hs,H,pl,q,p,volfrac);
    %%%% The parameters low, upp, raa0 and raa are calculated:
    [low,upp,raa0,raa] = ...
    MMA_asymp(outeriter,n,x,xold1,xold2,xmin,xmax,low,upp, ...
    raa0,raa,raa0eps,raaeps,df0dx,dfdx);
    [xmma,ymma,zmma,lam,xsi,eta,mu,zet,s,f0app,fapp] = ...
    MMA_gcmmasub(m,n,outeriter,epsimin,x,xmin,xmax,low,upp, ...
    raa0,raa,f0val,df0dx,fval,dfdx,a0,a,c,d,gf);
    xold2 = xold1;
    xold1 = x;
    x  = xmma;
    x_his(:,outit)=xmma;

    % PRINT RESULTS
    fprintf(' It.:%5i   P-norm Stress.:%11.4f   Vol.:%7.3f   MISES(max).:%11.4f \n',outit,pnorm, ...
        mean(x(:)),max(MISES));

    % Plot Density
    plot_x = -ceil(max(0,x(:)-plotCut));
    figure(f1); 
    patch('Faces',ELEM,'Vertices',NODE,'FaceVertexCData',-x,'FaceColor','flat','LineStyle','none'); 
    figure(f2); 
    patch('Faces',ELEM,'Vertices',NODE,'FaceVertexCData',MISES,'FaceColor','flat','LineStyle','none'); caxis([min(MISES) max(MISES)]);
    figure(f3); 
    patch('Faces',ELEM,'Vertices',NODE,'FaceVertexCData',plot_x,'FaceColor','flat','LineStyle','none'); 
    
    %%%% The residual vector of the KKT conditions is calculated:
    % [residu,kktnorm,residumax] = ...
    % MMA_kktcheck(m,n,xmma,ymma,zmma,lam,xsi,eta,mu,zet,s, ...
    % xmin,xmax,df0dx,fval,dfdx,a0,a,c,d);
%    outvector1 = [outeriter innerit x'];
%    outvector2 = [f0val fval'];

end
saveFunction(saveFileName,x,f1,f3,plotCut)
