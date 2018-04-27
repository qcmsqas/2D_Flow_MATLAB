function cal_XYAB(Nx,Ny)
    global Xf Yf Xv Yu XBu XBv YBu YBv XAu YAv;
    XBu = zeros(Nx+1,1);
    XAu = zeros(Nx+1,1);
    XBv = zeros(Nx+2,1);
    YBu = zeros(Ny+2,1);
    YBv = zeros(Ny+1,1);
    YAv = zeros(Ny+1,1);
    for i=2:Nx
        XBu(i)=1/(Xf(i+1)-Xf(i-1));
    end
    for i=2:Nx+1
        XBv(i)=1/(Xv(i+1)-Xv(i-1));
    end
    for i=2:Ny+1
        YBu(i)=1/(Yu(i+1)-Yu(i-1));
    end
    for i=2:Ny
        YBv(i)=1/(Yf(i+1)-Yf(i-1));
    end
    for i=2:Nx
        XAu(i)=1/(Xv(i+1)-Xv(i));
    end
    for i=2:Ny
        YAv(i)=1/(Yu(i+1)-Yu(i));
    end
    
end