function solveFlow(Nx,Ny,ux0,rmu,tmax,dt)
    %global u v p us vs Su Sv Mu Mv Mpu Mpv uus vvs XAu YAv XBu XBv YBu YBv fi;
    global fi p;
    cal_uvp_matrix(Nx,Ny,rmu,dt);
    cal_XYAB(Nx,Ny);
    oldfi = p(:);
    for it=0:dt:tmax
        cal_uv_source(Nx,Ny,ux0,dt);
        cal_usvs();
        cal_fi();
        cal_uv(Nx,Ny);
        err = max(max(abs(oldfi-fi)))/max(max(abs(fi)));
        disp(['calculating ' num2str(it) ' s, convergence = ' num2str(err)]);
        
        if err<1e-9
            break;
        end
        oldfi=fi;
    end
end