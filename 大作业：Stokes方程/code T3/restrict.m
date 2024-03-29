function [Fu,Fv,Fp]=restrict(N,fu,fv,fp)
    %fu:N*(N-1) fv:(N-1)*N  fp:N*N  
    Fu=1/4*fu(1:2:N-1,2:2:N-2)+1/4*fu(2:2:N,2:2:N-2);
    Fu=Fu+1/8*fu(1:2:N-1,1:2:N-3)+1/8*fu(1:2:N-1,3:2:N-1)+1/8*fu(2:2:N,1:2:N-3)+1/8*fu(2:2:N,3:2:N-1);
    Fv=1/4*fv(2:2:N-2,1:2:N-1)+1/4*fv(2:2:N-2,2:2:N);
    Fv=Fv+1/8*fv(1:2:N-3,1:2:N-1)+1/8*fv(3:2:N-1,1:2:N-1)+1/8*fv(1:2:N-3,2:2:N)+1/8*fv(3:2:N-1,2:2:N);
    Fp=1/4*fp(1:2:N-1,1:2:N-1)+1/4*fp(1:2:N-1,2:2:N)+1/4*fp(2:2:N,1:2:N-1)+1/4*fp(2:2:N,2:2:N);
end