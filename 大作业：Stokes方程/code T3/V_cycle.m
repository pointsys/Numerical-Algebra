function [u,v]=V_cycle(N,f,g,nu1,nu2,bottom)
    u=zeros(N,N+1);
    v=zeros(N+1,N);
    p=zeros(N,N);
    R=1;
    r0=sqrt(norm(f,'fro')^2+norm(g,'fro')^2);
    ru1=f;rv1=g;rp1=zeros(N,N);

    cnt=0;

    while R>=10^(-3)
        if cnt == 0
            m=1;
        else
            m=0;
        end
        [u1,v1]=vcycle(ru1,rv1,rp1,N,m,bottom,nu1,nu2);
        u=u+u1;v=v+v1;
        [u,v]=SymGS(N,u,v,p,f,g,nu2);
        ru1=ru(N,u,p,f);
        rv1=rv(N,v,p,g);
        r=sqrt(norm(ru1,'fro')^2+norm(rv1,'fro')^2);
        R=r/r0;
        cnt=cnt+1;
    end
end