eN=[];%用于存储误差
Cnt=[];
nu1=2;
nu2=2;
bottom=2;
tau=0.00001;
alpha=1;

for n = 6:11
    N=2^n;
    h=1/N;
    [f,g,u_true,v_true,p_true]=discrete(N);
    
    u=zeros(N,N+1);
    v=zeros(N+1,N);
    p=zeros(N,N);
    ru1=ru(N,u,p,f);
    rv1=rv(N,v,p,g);
    rp1=rp(N,u,v,zeros(N,N));
    
    R=1;
    r0=sqrt(norm(ru1,'fro')^2+norm(rv1,'fro')^2+norm(rp1,'fro')^2);
    
    cnt=0;
    
    tic;
    while R>=10^(-8)
        er=tau*norm(BT(N,u,v),'fro');
        [u,v]=PCG(N,f,g,p,er,1e-8,2*N*(N-1),bottom,nu1,nu2);

        p=p+alpha*BT(N,u,v);

        r=sqrt(norm(ru(N,u,p,f),'fro')^2 ...
            +norm(rv(N,v,p,g),'fro')^2 ...
            +norm(rp(N,u,v,zeros(N,N)),'fro')^2);
        R=r/r0;
        cnt=cnt+1;

    end

    toc;

    eu=u-u_true;
    ev=v-v_true;
    e1=sqrt(norm(eu,'fro')^2+norm(ev,'fro')^2)*h;
    eN=[eN,e1];
    Cnt=[Cnt,cnt];
end