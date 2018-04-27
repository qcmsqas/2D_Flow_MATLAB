function setMesh(xmin,xmax,ymin,ymax,Nx,Ny)
    global p u v us vs Xc Xf Yc Yf Xv Yu;
    p=zeros(Nx,Ny);
    u=zeros(Nx+1,Ny+2);
    v=zeros(Nx+2,Ny+1);
    us=u(:);
    vs=v(:);
    dx=(xmax-xmin)/Nx;
    Xf=xmin:dx:xmax+dx*0.1;
    Xc=xmin+dx*0.5:dx:xmax;
    dy=(ymax-ymin)/Ny;
    Yf=ymin:dy:ymax+dx*0.1;
    Yc=ymin+dy*0.5:dy:ymax;
    Xv=xmin-dx*0.5:dx:xmax+dx*0.6;
    Yu=ymin-dy*0.5:dy:ymax+dx*0.6;
end