nu1=3;
nu2=3;
bottom=4;

eN=[];%用于存储误差
Cnt=[];%用于存储V-cycle次数

for n = 6:11
    N=2^n;
    h=1/N;
    [f,g,u_true,v_true,p_true]=discrete(N);

    tic;
    
    u=zeros(N,N+1);
    v=zeros(N+1,N);
    p=zeros(N,N);
    R=1;
    r0=sqrt(norm(f,'fro')^2+norm(g,'fro')^2);
    ru1=f;rv1=g;rp1=zeros(N,N);

    cnt=0;

    while R>=10^(-8)
        if cnt == 0
            m=1;
        else
            m=0;
        end
        [u1,v1,p1]=vcycle(ru1,rv1,rp1,N,m,bottom,nu1,nu2);
        u=u+u1;v=v+v1;p=p+p1;
        [u,v,p]=DGS(N,u,v,p,f,g,zeros(N,N),nu2);
        ru1=ru(N,u,p,f);
        rv1=rv(N,v,p,g);
        rp1=rp(N,u,v,zeros(N,N));
        r=sqrt(norm(ru1,'fro')^2+norm(rv1,'fro')^2+norm(rp1,'fro')^2);
        R=r/r0;
        cnt=cnt+1;
    end
    toc;

    Cnt=[Cnt,cnt];

    eu=u-u_true;
    ev=v-v_true;
    e1=sqrt(norm(eu,'fro')^2+norm(ev,'fro')^2)*h;

    eN=[eN,e1];
end