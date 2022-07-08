% Main_function(0.3,3,0.5,10,0.3)   pl:penal, q:stress relaxation, p:p-norm
function Main_function(rmin,pl,q,p,volfrac)
addpath('FE'); addpath('MMA'); addpath('data');
[NODE,ELEM] = inp_('main.inp');
[Hs,H]=Prepare_filter(rmin,NODE,ELEM);
nele = length(ELEM);
x=volfrac*ones(nele,1);
outeriter = 0;
%%%%%%%%%%%%%%%%%% M M A Zone %%%%%%%%%%%%%%%%%%
m =1;
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
c       = [1e6]';
d       = [0]';
a0      = 1;
a       = [0]';
raa0    = 0.0001;
raa     = 0.0001;
raa0eps = 0.0000001;
raaeps  = 0.0000001;
maxoutit  = 120;
kkttol  = 0;
nele = length(ELEM);
x_his=zeros(nele,maxoutit);
kktnorm = kkttol+1;
outit = 0;
while  outit < maxoutit
    outit   = outit+1;
    if outit <= 15, gf = 0.2; else gf = min(0.5,1.01*gf); end
    outeriter = outeriter+1;
    [f0val,df0dx,fval,dfdx,MISES_MAX]=Stress_minimize(NODE,ELEM,x,Hs,H,pl,q,p,volfrac);
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

    % PRINT RESULTS
    fprintf(' It.:%5i   P-norm Stress.:%11.4f   Vol.:%7.3f   MISES(max).:%11.4f \n',outit,f0val, ...
        mean(x(:)),MISES_MAX);
    %%%% The residual vector of the KKT conditions is calculated:
    [residu,kktnorm,residumax] = ...
    MMA_kktcheck(m,n,xmma,ymma,zmma,lam,xsi,eta,mu,zet,s, ...
    xmin,xmax,df0dx,fval,dfdx,a0,a,c,d);
%    outvector1 = [outeriter innerit x'];
%    outvector2 = [f0val fval'];
    x_his(:,outit)=xmma;
end
saveX=fliplr(reshape(x,[100,100]))';
save('1.mat','saveX')
