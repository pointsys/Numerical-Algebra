eN=[];%用于存储误差
for n = 6:9
    N=2^n;
    h=1/N;
    [f,g,u_true,v_true,p_true]=discrete(N);
    
    eu0=zeros(N,N+1);
    ev0=zeros(N+1,N);
    ep0=zeros(N,N);
    ru1=ru(N,eu0,ep0,f);
    rv1=rv(N,ev0,ep0,g);
    rp1=rp(N,eu0,ev0,zeros(N,N));
    
    R=1;
    r0=sqrt(norm(ru1,'fro')^2+norm(rv1,'fro')^2+norm(rp1,'fro')^2);
    
    p=zeros(N,N);

    A=FacA(N);%对A做带状Gauss消去
    while R>=10^(-8)

        [u,v]=divA(N,A,p,f,g);

        p=p+BT(N,u,v);

        r=sqrt(norm(ru(N,u,p,f),'fro')^2 ...
            +norm(rv(N,v,p,g),'fro')^2 ...
            +norm(rp(N,u,v,zeros(N,N)),'fro')^2);
        R=r/r0;

    end

    eu=u-u_true;
    ev=v-v_true;
    e1=sqrt(norm(eu,'fro')^2+norm(ev,'fro')^2)*h;
    eN=[eN,e1];
end