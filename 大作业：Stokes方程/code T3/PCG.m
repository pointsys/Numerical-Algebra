function [u,v]=PCG(N,f,g,p,er,epsilon,kmax,bottom,nu1,nu2)
    k=0;
    f=f-Bup(N,p);g=g-Bvp(N,p);
    ru1=f;rv1=g;
    u=zeros(N,N+1);
    v=zeros(N+1,N);
    b=sqrt(norm(ru1,'fro')^2+norm(rv1,'fro')^2);
    while sqrt(norm(ru1,'fro')^2+norm(rv1,'fro')^2)>max(epsilon*b,er) && k<kmax
        k=k+1;
        [u_z,v_z]=V_cycle(N,ru1,rv1,nu1,nu2,bottom);
        if k == 1
            p_u=u_z;p_v=v_z;
            rho=sum(ru1.*u_z(:,2:N),'all')+sum(rv1.*v_z(2:N,:),'all');
        else
            rho_old=rho;
            rho=sum(ru1.*u_z(:,2:N),'all')+sum(rv1.*v_z(2:N,:),'all');
            beta=rho/rho_old;
            p_u=u_z+beta*p_u;p_v=v_z+beta*p_v;
        end
        u_w=Au(N,p_u,zeros(N,N));v_w=Av(N,p_v,zeros(N,N));
        alpha=rho/(sum(p_u(:,2:N).*u_w,'all')+sum(p_v(2:N,:).*v_w,"all"));
        u=u+alpha*p_u;v=v+alpha*p_v;
        ru1=ru1-alpha*u_w;rv1=rv1-alpha*v_w;
        
    end
end