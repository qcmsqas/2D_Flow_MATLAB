function cal_uv(Nx,Ny)
    global XAu YAv us vs u v fi;
    nnx=Nx+1;
    for i=2:Nx
        for jj=2:Ny+1
            j=jj-1;
            u(i,jj)=us(i+j*nnx)-(fi(i+(j-1)*Nx)-fi(i-1+(j-1)*Nx))*XAu(i);
        end
    end
    nnx=Nx+2;
    for i=2:Nx+1
        for jj=2:Ny
            j=jj-1;
            v(i,jj)=vs(i+j*nnx)-(fi(i-1+j*Nx)-fi(i-1+(j-1)*Nx))*YAv(jj);
        end
    end
end