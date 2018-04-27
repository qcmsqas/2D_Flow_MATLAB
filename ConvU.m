function [UU]=ConvU(U)
    N=size(U,1);
    M=size(U,3);
    UU=zeros(N*M,N*M);
    for i=1:M-1
    UU(1+(i-1)*N:N+(i-1)*N,1+(i)*N:N+(i)*N)=U(:,:,i);
    end


end