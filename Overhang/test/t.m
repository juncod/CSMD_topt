if nSens
    dfxi=varargin; dfx=varargin; 
    lambda = zeros(nSens, nelx,nelz);
    for i=1:nely-1
        dsmindx  = .5*(1-(x(i,:,:)-Xi(i,:,:))./sq(i,:,:));
        dsmindXi = 1-dsmindx; 
        cbr = zeros(1,nelx+2,nelz+2);
        cbr(1,2:nelx+1, 2:nelz+1) = xi(i+1,:,:)+SHIFT;
        dmx = zeros(Ns, nelx, nelz);
        for j = 1:Ns
            if j <= 3
                dmx(j,:,:) = (P/Q)*keep(i,:,:).^(1/Q-1).*cbr(1,(1:nelx)+(j-1),(2:nelz+1)).^(P-1);
            elseif j == 4
                dmx(j,:,:) = (P/Q)*keep(i,:,:).^(1/Q-1).*cbr(1,(2:nelx+1),(1:nelz)).^(P-1);
            elseif j == 5
                dmx(j,:,:) = (P/Q)*keep(i,:,:).^(1/Q-1).*cbr(1,(2:nelx+1),(1:nelz)+2).^(P-1);
            end
        end
        for k=1:nSens
            dfx{k}(i,:,:) = dsmindx.*(dfxi{k}(i,:,:)+lambda(k,:,:));
            preLambda = zeros(1,nelx+2,nelz+2);
            preLambda(1,2:nelx+1,2:nelz+1) = (dfxi{k}(i,:,:)+lambda(k,:,:)).*dsmindXi;
            for nz=1:nelz
                for nx = 1:nelx
                    lambda(k, nx, nz) = sum([preLambda(1, nx-1+1,nz+1);...
                                            preLambda(1, nx+1,nz+1);...
                                            preLambda(1, nx+1+1,nz+1);...
                                            preLambda(1, nx+1,nz-1+1);...
                                            preLambda(1, nx+1,nz+1+1)].*...
                                            dmx(2,nx,nz));
                end
            end
        end
    end
