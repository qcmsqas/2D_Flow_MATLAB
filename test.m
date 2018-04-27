 err = u(2:end,2:end-1)-u(1:end-1,2:end-1)+v(2:end-1,2:end)-v(2:end-1,1:end-1)
pout=getP(Nx,Ny);
pout=pout*dt/ro;
ppout=zeros(Nx+2,Ny+2);
ppout(2:end-1,2:end-1)=pout;
    ppout(:,1)=ppout(:,2);
    ppout(:,end)=ppout(:,end-1);
    ppout(1,:)=ppout(2,:);
    ppout(end,:)=ppout(end-1,:);
    
Sp = MMp*fi;
ssp = Sp;
for i=1:Nx
    for jj=1:Ny
        j=jj-1;
        ssp(i+j*Nx)=us(i+1+(j+1)*(Nx+1))-us(i+(j+1)*(Nx+1))+vs(i+1+(j+1)*(Nx+2))-vs(i+1+j*(Nx+2));
    end
end
dx=Xc(2)-Xc(1);
dy=dx;
ssp=ssp/dx;