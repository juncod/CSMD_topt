if nSens
    dfxi=varargin; dfx=varargin; 
    lambda=zeros(nSens,nelx,nelz); 
    for i=1:nely-1
        dsmindx  = .5*(1-(x(i,:,:)-Xi(i,:,:))./sq(i,:,:));
        dsmindXi = 1-dsmindx; 
        cbr = padarray(squeeze(xi(i+1,:,:)),[1 1],0,'both');
        dmx = zeros(Ns,nelx,nelz);
        for j=1:Ns
            dmx(j,:,:) = (P/Q)*keep(i,:,:).^(1/Q-1).*cbr((1:nelx)+(j-1)).^(P-1);
        end        
        for k=1:nSens
            dfx{k}(i,:,:) = dsmindx.*(dfxi{k}(i,:,:)+lambda(k,:,:));
            for lx=1:nelx
                for lz=1:nelz
                    dsmaxdxi = zeros(nelx+2,nelz+2);
                    dsmaxdxi(lx,lz+1) = dmx(1,lx,lz);
                    dsmaxdxi(lx+1,lz+1) = dmx(2,lx,lz);
                    dsmaxdxi(lx+2,lz+1) = dmx(3,lx,lz);
                    dsmaxdxi(lx+1,lz) = dmx(4,lx,lz);
                    dsmaxdxi(lx+1,lz+2) = dmx(5,lx,lz);
                    dsmaxdxi(end,:)=[];dsmaxdxi(:,end)=[]; dsmaxdxi(1,:)=[]; dsmaxdxi(:,1)=[];
                    dsmaxdxi=reshape(dsmaxdxi,1,nelx,nelz);

                    lambda(k,lx,lz)= sum(sum(((dfxi{k}(i,:,:)+lambda(k,:,:)).*dsmindXi).*dsmaxdxi));
                    end
                end
            end
    end
            