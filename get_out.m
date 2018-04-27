function [pout,uout,vout,uvout]=get_out(Nx,Ny,ro,dt)
    global fi u v;
    
    pout=zeros(Nx,Ny);
    uout=pout;
    vout=pout;
    uvout=pout;
    for i=1:Nx
        for j=1:Ny
            pout(i,j)=fi(i+(j-1)*Nx);
            uout(i,j)=(u(i+1,j+1)+u(i,j+1))/2;
            vout(i,j)=(v(i+1,j+1)+v(i+1,j))/2;
            uvout(i,j)=(uout(i,j)^2+vout(i,j)^2)^0.5;
        end
    end
    pout=pout*ro/dt;
end