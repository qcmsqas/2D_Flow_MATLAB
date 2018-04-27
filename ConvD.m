function [DD]=ConvD(D)
    N=size(D,1);
    M=size(D,3);
    DD=zeros(N*M,N*M);
    for i=1:M
    DD(1+(i-1)*N:N+(i-1)*N,1+(i-1)*N:N+(i-1)*N)=D(:,:,i);
    end


end