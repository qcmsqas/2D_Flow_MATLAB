function [LL]=ConvL(L)
    N=size(L,1);
    M=size(L,3);
    LL=zeros(N*M,N*M);
    for i=1:M-1
    LL(1+(i)*N:N+(i)*N,1+(i-1)*N:N+(i-1)*N)=L(:,:,i+1);
    end


end