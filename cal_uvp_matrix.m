function cal_uvp_matrix(Nx,Ny,rmu,dt)
    q=1/dt;
    global Xc Yc Xv Yu Xf Yf;
    %获得u更新us的泊松方程的方程矩阵。
    global Mu;
    Mu=zeros((Nx+1)*(Ny+2));
    nnx=Nx+1;
    %%考虑内部点
    for i=2:Nx
        for jj=2:Ny+1
            j=jj-1;
            Mu(i+j*nnx,i+j*nnx)=Mu(i+j*nnx,i+j*nnx)+(q ...
                              +rmu/(Xf(i)-Xf(i-1))/(Xc(i)-Xc(i-1)) ...
                              +rmu/(Xf(i+1)-Xf(i))/(Xc(i)-Xc(i-1)) ...
                              +rmu/(Yu(jj)-Yu(jj-1))/(Yf(jj)-Yf(jj-1)) ...
                              +rmu/(Yu(jj+1)-Yu(jj))/(Yf(jj)-Yf(jj-1)) ...
                              );
            Mu(i+j*nnx,(i-1)+j*nnx)=Mu(i+j*nnx,(i-1)+j*nnx)-rmu/(Xf(i)-Xf(i-1))/(Xc(i)-Xc(i-1));
            Mu(i+j*nnx,(i+1)+j*nnx)=Mu(i+j*nnx,(i+1)+j*nnx)-rmu/(Xf(i+1)-Xf(i))/(Xc(i)-Xc(i-1));
            Mu(i+j*nnx,i+(j-1)*nnx)=Mu(i+j*nnx,i+(j-1)*nnx)-rmu/(Yu(jj)-Yu(jj-1))/(Yf(jj)-Yf(jj-1));
            Mu(i+j*nnx,i+(j+1)*nnx)=Mu(i+j*nnx,i+(j+1)*nnx)-rmu/(Yu(jj+1)-Yu(jj))/(Yf(jj)-Yf(jj-1));
        end
    end
    %%考虑上下边界
    for i=1:Nx+1
        %%%下边界
        Mu(i,i)=Mu(i,i)+1;
        Mu(i,i+nnx)=Mu(i,i+nnx)+1;
        %%%上边界
        Mu(i+(Ny+1)*nnx,i+(Ny+1)*nnx)=Mu(i+(Ny+1)*nnx,i+(Ny+1)*nnx)+1;
        Mu(i+(Ny+1)*nnx,i+(Ny)*nnx)=Mu(i+(Ny+1)*nnx,i+(Ny)*nnx)+1;
    end
    %%考虑左右边界
    for jj=2:Ny+1
        j=jj-1;
        Mu(1+j*nnx,1+j*nnx)=Mu(1+j*nnx,1+j*nnx)+1;
        Mu(Nx+1+j*nnx,Nx+1+j*nnx)=Mu(Nx+1+j*nnx,Nx+1+j*nnx)+1;
    end
    Mu=Mu^(-1);
    
    
    %获得v更新vs的泊松方程的方程矩阵
    global Mv;
    Mv=zeros((Nx+2)*(Ny+1));
    nnx=Nx+2;
    %%考虑内部点
    for i=2:Nx+1
        for jj=2:Ny
            j=jj-1;
            Mv(i+j*nnx,i+j*nnx)=Mv(i+j*nnx,i+j*nnx)+(q ...
                       +rmu/(Xv(i)-Xv(i-1))/(Xf(i)-Xf(i-1)) ...
                       +rmu/(Xv(i+1)-Xv(i))/(Xf(i)-Xf(i-1)) ...
                       +rmu/(Yf(jj)-Yf(jj-1))/(Yc(jj)-Yc(jj-1)) ...
                       +rmu/(Yf(jj+1)-Yf(jj))/(Yc(jj)-Yc(jj-1)) ...
                       );
            Mv(i+j*nnx,(i-1)+j*nnx)=Mv(i+j*nnx,(i-1)+j*nnx)-rmu/(Xv(i)-Xv(i-1))/(Xf(i)-Xf(i-1));
            Mv(i+j*nnx,(i+1)+j*nnx)=Mv(i+j*nnx,(i+1)+j*nnx)-rmu/(Xv(i+1)-Xv(i))/(Xf(i)-Xf(i-1));
            Mv(i+j*nnx,i+(j-1)*nnx)=Mv(i+j*nnx,i+(j-1)*nnx)-rmu/(Yf(jj)-Yf(jj-1))/(Yc(jj)-Yc(jj-1));
            Mv(i+j*nnx,i+(j+1)*nnx)=Mv(i+j*nnx,i+(j+1)*nnx)-rmu/(Yf(jj+1)-Yf(jj))/(Yc(jj)-Yc(jj-1));
        end
    end
    %%考虑左右边界
    for jj=1:Ny+1
        j=jj-1;
        %%%左边界
        Mv(1+j*nnx,1+j*nnx)=Mv(1+j*nnx,1+j*nnx)+1;
        Mv(1+j*nnx,2+j*nnx)=Mv(1+j*nnx,2+j*nnx)+1;
        %%%%右边界
        Mv(Nx+2+j*nnx,Nx+2+j*nnx)=Mv(Nx+2+j*nnx,Nx+2+j*nnx)+1;
        Mv(Nx+2+j*nnx,Nx+1+j*nnx)=Mv(Nx+2+j*nnx,Nx+1+j*nnx)+1;
    end
    %%考虑上下边界
    for i=2:Nx+1
        %%%下边界
        Mv(i,i)=Mv(i,i)+1;
        %%%上边界
        Mv(i+Ny*nnx,i+Ny*nnx)=Mv(i+Ny*nnx,i+Ny*nnx)+1;
    end
    Mv=Mv^(-1);
    
    %计算p的泊松方程的方程矩阵
    global Mpu Mpv;
    Mp=zeros(Nx*Ny);
    %%考虑内部点
    for i=1:Nx
        for jj=1:Ny
            j=jj-1;
            %disp([num2str(i) ' ' num2str(jj)]);
            %处理四种边界
            ip=1;iq=1;jp=1;jq=1;
            %必须有一个狄利克雷边界，否则会不稳定
            if i==1&&jj==1
                Mp(i+j*Nx,i+j*Nx)=Mp(i+j*Nx,i+j*Nx)+4/(Xv(2)-Xv(1))/(Yu(2)-Yu(1));
                continue;
            end
            if i==1
                ip=0;
            end
            if i==Nx
                iq=0;
            end
            if jj==1
                jp=0;
            end
            if jj==Ny
                jq=0;
            end
            Mp(i+j*Nx,i+j*Nx)=Mp(i+j*Nx,i+j*Nx)-1/(Xv(i+1)-Xv(i))/(Xf(i+1)-Xf(i)) ...
                                               -1/(Xv(i+2)-Xv(i+1))/(Xf(i+1)-Xf(i)) ...
                                               -1/(Yu(jj+1)-Yu(jj))/(Yf(jj+1)-Yf(jj)) ...
                                               -1/(Yu(jj+2)-Yu(jj+1))/(Yf(jj+1)-Yf(jj));
            Mp(i+j*Nx,(i-ip)+j*Nx)=Mp(i+j*Nx,(i-ip)+j*Nx)+1/(Xv(i+1)-Xv(i))/(Xf(i+1)-Xf(i));
            Mp(i+j*Nx,(i+iq)+j*Nx)=Mp(i+j*Nx,(i+iq)+j*Nx)+1/(Xv(i+2)-Xv(i+1))/(Xf(i+1)-Xf(i));
            Mp(i+j*Nx,i+(j-jp)*Nx)=Mp(i+j*Nx,i+(j-jp)*Nx)+1/(Yu(jj+1)-Yu(jj))/(Yf(jj+1)-Yf(jj));
            Mp(i+j*Nx,i+(j+jq)*Nx)=Mp(i+j*Nx,i+(j+jq)*Nx)+1/(Yu(jj+2)-Yu(jj+1))/(Yf(jj+1)-Yf(jj));
        end
    end
    
    Mp=Mp^(-1);

    MMpu=zeros(Nx*Ny,(Nx+1)*(Ny+2));
    MMpv=zeros(Nx*Ny,(Nx+2)*(Ny+1));
    for i=1:Nx
        for jj=1:Ny
            %为了方程稳定性，必须有一个狄利克雷边条
            if i==1&&jj==1
                continue;
            end
            j=jj-1;
            MMpu(i+j*Nx,i+(j+1)*(Nx+1))=-1/(Xf(i+1)-Xf(i));
            MMpu(i+j*Nx,i+1+(j+1)*(Nx+1))=1/(Xf(i+1)-Xf(i));
            MMpv(i+j*Nx,i+1+j*(Nx+2))=-1/(Yf(jj+1)-Yf(jj));
            MMpv(i+j*Nx,i+1+(j+1)*(Nx+2))=1/(Yf(jj+1)-Yf(jj));
        end
    end
    Mpu=Mp*MMpu;
    Mpv=Mp*MMpv;
    
end