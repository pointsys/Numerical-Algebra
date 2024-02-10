function [u,v,p]=DGS(N,u,v,p,f,g,d,nu)
%进行nu步DGS迭代

    h=1/N;
    for k = 1:nu
        
        %用Gauss-Seidel迭代更新速度分量
        ru1=ru(N,u,p,f);

        for i = 1:N-2
            ru1(1,i)=ru1(1,i)*h^2/3;
            ru1(1,i+1)=ru1(1,i+1)+(1/h^2)*ru1(1,i);
            ru1(2,i)=ru1(2,i)+(1/h^2)*ru1(1,i);
            for j = 2:N-1
                ru1(j,i)=ru1(j,i)*h^2/4;
                ru1(j,i+1)=ru1(j,i+1)+(1/h^2)*ru1(j,i);
                ru1(j+1,i)=ru1(j+1,i)+(1/h^2)*ru1(j,i);
            end
            ru1(N,i)=ru1(N,i)*h^2/3;
            ru1(N,i+1)=ru1(N,i+1)+(1/h^2)*ru1(N,i);
        end

        ru1(1,N-1)=ru1(1,N-1)*h^2/3;
        ru1(2,N-1)=ru1(2,N-1)+(1/h^2)*ru1(1,N-1);

        for j = 2:N-1
            ru1(j,N-1)=ru1(j,N-1)*h^2/4;
            ru1(j+1,N-1)=ru1(j+1,N-1)+(1/h^2)*ru1(j,N-1);
        end
        ru1(N,N-1)=ru1(N,N-1)*h^2/3;
        
        u(:,2:N)=u(:,2:N)+ru1;

        rv1=rv(N,v,p,g);

        for i = 1:N-2
            rv1(i,1)=rv1(i,1)*h^2/3;
            rv1(i+1,1)=rv1(i+1,1)+(1/h^2)*rv1(i,1);
            rv1(i,2)=rv1(i,2)+(1/h^2)*rv1(i,1);
            for j = 2:N-1
                rv1(i,j)=rv1(i,j)*h^2/4;
                rv1(i+1,j)=rv1(i+1,j)+(1/h^2)*rv1(i,j);
                rv1(i,j+1)=rv1(i,j+1)+(1/h^2)*rv1(i,j);
            end
            rv1(i,N)=rv1(i,N)*h^2/3;
            rv1(i+1,N)=rv1(i+1,N)+(1/h^2)*rv1(i,N);
        end
        
        rv1(N-1,1)=rv1(N-1,1)*h^2/3;
        rv1(N-1,2)=rv1(N-1,2)+(1/h^2)*rv1(N-1,1);
        for j = 2:N-1
            rv1(N-1,j)=rv1(N-1,j)*h^2/4;
            rv1(N-1,j+1)=rv1(N-1,j+1)+(1/h^2)*rv1(N-1,j);
        end
        rv1(N-1,N)=rv1(N-1,N)*h^2/3;
        
        v(2:N,:)=v(2:N,:)+rv1;
        
        for i = 2:N-1
            for j = 2:N-1
                rp1=-1/h*(u(j,i+1)-u(j,i)+v(j+1,i)-v(j,i))-d(j,i);
                delta=rp1*h/4;
                %更新内部单元速度
                u(j,i)=u(j,i)-delta;
                u(j,i+1)=u(j,i+1)+delta;
                v(j,i)=v(j,i)-delta;
                v(j+1,i)=v(j+1,i)+delta;
                %更新内部单元压力
                p(j,i)=p(j,i)+rp1;
                p(j+1,i)=p(j+1,i)-rp1/4;
                p(j-1,i)=p(j-1,i)-rp1/4;
                p(j,i+1)=p(j,i+1)-rp1/4;
                p(j,i-1)=p(j,i-1)-rp1/4;
            end
        end
        
        
        for i = 2:N-1
            rp1=-1/h*(u(N,i+1)-u(N,i)+v(N+1,i)-v(N,i))-d(j,i);
            delta=rp1*h/3;
            %更新边界单元j=N速度
            u(N,i)=u(N,i)-delta;
            u(N,i+1)=u(N,i+1)+delta;
            v(N,i)=v(N,i)-delta;
            %更新边界单元j=N压力
            p(N,i)=p(N,i)+rp1;
            p(N,i+1)=p(N,i+1)-rp1/3;
            p(N,i-1)=p(N,i-1)-rp1/3;
            p(N-1,i)=p(N-1,i)-rp1/3;
        end
        
        
        for j = 2:N-1
            rp1=-1/h*(u(j,N+1)-u(j,N)+v(j+1,N)-v(j,N))-d(j,N);
            delta=rp1*h/3;
            %更新边界单元i=N速度
            v(j,N)=v(j,N)-delta;
            v(j+1,N)=v(j+1,N)+delta;
            u(j,N)=u(j,N)-delta;
            %更新边界单元i=N压力
            p(j,N)=p(j,N)+rp1;
            p(j+1,N)=p(j+1,N)-rp1/3;
            p(j-1,N)=p(j-1,N)-rp1/3;
            p(j,N-1)=p(j,N-1)-rp1/3;
        end
        
        
        for j = 2:N-1
            rp1=-1/h*(u(j,2)-u(j,1)+v(j+1,1)-v(j,1))-d(j,1);
            delta=rp1*h/3;
            %更新边界单元i=1速度
            v(j,1)=v(j,1)-delta;
            v(j+1,1)=v(j+1,1)+delta;
            u(j,2)=u(j,2)+delta;
            %更新边界单元i=1压力
            p(j,1)=p(j,1)+rp1;
            p(j+1,1)=p(j+1,1)-rp1/3;
            p(j-1,1)=p(j-1,1)-rp1/3;
            p(j,2)=p(j,2)-rp1/3;
        end
        
        
        for i = 2:N-1
            rp1=-1/h*(u(1,i+1)-u(1,i)+v(2,i)-v(1,i))-d(1,i);
            delta=rp1*h/3;
            %更新边界单元j=1速度
            u(1,i)=u(1,i)-delta;
            u(1,i+1)=u(1,i+1)+delta;
            v(2,i)=v(2,i)+delta;
            %更新边界单元j=1压力
            p(1,i)=p(1,i)+rp1;
            p(1,i+1)=p(1,i+1)-rp1/3;
            p(1,i-1)=p(1,i-1)-rp1/3;
            p(2,i)=p(2,i)-rp1/3;
        end
        
        rp1=-1/h*(u(1,2)-u(1,1)+v(2,1)-v(1,1))-d(1,1);
        %更新顶点单元(1,1)速度
        delta=rp1*h/2;
        u(1,2)=u(1,2)+delta;
        v(2,1)=v(2,1)+delta;
        %更新顶点单元(1,1)压力
        p(1,1)=p(1,1)+rp1;
        p(2,1)=p(2,1)-rp1/2;
        p(1,2)=p(1,2)-rp1/2;

        rp1=-1/h*(u(N,2)-u(N,1)+v(N+1,1)-v(N,1))-d(N,1);
        %更新顶点单元(1,N)速度
        delta=rp1*h/2;
        u(N,2)=u(N,2)+delta;
        v(N,1)=v(N,1)-delta;
        %更新顶点单元(1,N)压力
        p(N,1)=p(N,1)+rp1;
        p(N-1,1)=p(N-1,1)-rp1/2;
        p(N,2)=p(N,2)-rp1/2;

        rp1=-1/h*(u(1,N+1)-u(1,N)+v(2,N)-v(1,N))-d(1,N);
        %更新顶点单元(N,1)速度
        delta=rp1*h/2;
        v(2,N)=v(2,N)+delta;
        u(1,N)=u(1,N)-delta;
        %更新顶点单元(N,1)压力
        p(1,N)=p(1,N)+rp1;
        p(1,N-1)=p(1,N-1)-rp1/2;
        p(2,N)=p(2,N)-rp1/2;

        rp1=-1/h*(u(N,N+1)-u(N,N)+v(N+1,N)-v(N,N))-d(N,N);
        %更新顶点单元(N,N)速度
        delta=rp1*h/2;
        u(N,N)=u(N,N)-delta;
        v(N,N)=v(N,N)-delta;
        %更新顶点单元(N,N)压力
        p(N,N)=p(N,N)+rp1;
        p(N-1,N)=p(N-1,N)-rp1/2;
        p(N,N-1)=p(N,N-1)-rp1/2;

    end
end