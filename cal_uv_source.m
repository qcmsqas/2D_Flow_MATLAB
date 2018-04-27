function cal_uv_source(Nx,Ny,ux0,dt)
    global u v XBu XBv YBu YBv Su Sv;
    q=1/dt;
    %首先应用边界条件
    for i=1:Nx+1
        u(i,Ny+2)=2*ux0-u(i,Ny+1);
        u(i,1)=-u(i,2);
    end
    for j=1:Ny+1
        v(Nx+2,j)=-v(Nx+1,j);
        v(1,j)=-v(2,j);
    end
    %计算u方程源矩阵
    Su=zeros((Nx+1)*(Ny+2),1);
    nnx=Nx+1;
    %%计算内部情况
    for i=2:Nx
        for jj=2:Ny+1
            j=jj-1;
            Su(i+j*nnx)=-u(i,jj)*XBu(i)*(u(i+1,jj)-u(i-1,jj)) ... 
                       -0.25*(v(i,jj-1)+v(i,jj)+v(i+1,jj-1)+v(i+1,jj))*YBu(jj)*(u(i,jj+1)-u(i,jj-1)) ...
                       +q*u(i,jj);
        end
    end
    %%计算上下边界
    for i=1:Nx+1
        %%%下边界
        Su(i)=0;
        %%%上边界
        Su(i+(Ny+1)*nnx)=2*ux0;
    end
    %%计算左右边界
    for jj=2:Ny+1
        j=jj-1;
        Su(1+j*nnx)=0;
        Su(Nx+1+j*nnx)=0;
    end
    
    
    %计算v方程源矩阵
    Sv=zeros((Nx+2)*(Ny+1),1);
    nnx=Nx+2;
    %%计算内部情况
    for i=2:Nx+1
        for jj=2:Ny
            j=jj-1;
            Sv(i+j*nnx)=-0.25*(u(i-1,jj)+u(i-1,jj+1)+u(i,jj)+u(i,jj+1))*XBv(i)*(v(i+1,jj)-v(i-1,jj)) ... 
                       -v(i,jj)*YBv(jj)*(v(i,jj+1)-v(i,jj-1)) ...
                       +q*v(i,jj);
        end
    end
    %%计算上下边界
    for i=2:Nx+1
        %%%下边界
        Sv(i)=0;
        %%%上边界
        Sv(i+Ny*nnx)=0;
    end
    %%计算左右边界
    for jj=1:Ny+1
        j=jj-1;
        Sv(1+j*nnx)=0;
        Sv(Nx+2+j*nnx)=0;
    end
end