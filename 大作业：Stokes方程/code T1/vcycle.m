function [u,v,p]=vcycle(ru0,rv0,rp0,N,m,bottom,nu1,nu2)

    if N==bottom

        [u,v,p]=DGS(N,zeros(N,N+1),zeros(N+1,N),zeros(N,N),ru0,rv0,rp0,50);

    else
        if m == 1
            [u1,v1,p1]=DGS(N,zeros(N,N+1),zeros(N+1,N),zeros(N,N),ru0,rv0,rp0,nu1);%以0为初值，对残量方程进行迭代
            ru3=ru(N,u1,p1,ru0);rv3=rv(N,v1,p1,rv0);rp3=rp(N,u1,v1,rp0);%残量

            [ru1,rv1,rp1]=restrict(N,ru3,rv3,rp3);%限制残量
        else
            [ru1,rv1,rp1]=restrict(N,ru0,rv0,rp0);%限制残量
        end

        [eu2,ev2,ep2]=vcycle(ru1,rv1,rp1,N/2,1,bottom,nu1,nu2);%粗网格上误差，即细网格残量方程的近似解
        
        [u,v,p]=lift(N/2,eu2,ev2,ep2);%细网格上误差进行提升
        
        if m == 1
            u=u+u1;v=v+v1;p=p+p1;
        end
        
        [u,v,p]=DGS(N,u,v,p,ru0,rv0,rp0,nu2);%对残量方程再磨光
        

    end

end