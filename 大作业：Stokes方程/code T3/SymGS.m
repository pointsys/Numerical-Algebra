function [u,v]=SymGS(N,u,v,p,f,g,nu)%进行nu次对称迭代
    h=1/N;
    for k = 1:nu
        
        %第一步迭代
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
        
        %第二步迭代
        ru1=ru(N,u,p,f);

        for i = N-1:-1:2
            ru1(N,i)=ru1(N,i)*h^2/3;
            ru1(N,i-1)=ru1(N,i-1)+(1/h^2)*ru1(N,i);
            ru1(N-1,i)=ru1(N-1,i)+(1/h^2)*ru1(N,i);
            for j = N-1:-1:2
                ru1(j,i)=ru1(j,i)*h^2/4;
                ru1(j,i-1)=ru1(j,i-1)+(1/h^2)*ru1(j,i);
                ru1(j-1,i)=ru1(j-1,i)+(1/h^2)*ru1(j,i);
            end
            ru1(1,i)=ru1(1,i)*h^2/3;
            ru1(1,i-1)=ru1(1,i-1)+(1/h^2)*ru1(1,i);
        end

        ru1(N,1)=ru1(N,1)*h^2/3;
        ru1(N-1,1)=ru1(N-1,1)+(1/h^2)*ru1(N,1);

        for j = N-1:-1:2
            ru1(j,1)=ru1(j,1)*h^2/4;
            ru1(j-1,1)=ru1(j-1,1)+(1/h^2)*ru1(j,1);
        end
        ru1(1,1)=ru1(1,1)*h^2/3;
        
        u(:,2:N)=u(:,2:N)+ru1;

        rv1=rv(N,v,p,g);

        for j = N-1:-1:2
            rv1(j,N)=rv1(j,N)*h^2/3;
            rv1(j-1,N)=rv1(j-1,N)+(1/h^2)*rv1(j,N);
            rv1(j,N-1)=rv1(j,N-1)+(1/h^2)*rv1(j,N);
            for i = N-1:-1:2
                rv1(j,i)=rv1(j,i)*h^2/4;
                rv1(j,i-1)=rv1(j,i-1)+(1/h^2)*rv1(j,i);
                rv1(j-1,i)=rv1(j-1,i)+(1/h^2)*rv1(j,i);
            end
            rv1(j,1)=rv1(j,1)*h^2/3;
            rv1(j-1,1)=rv1(j-1,1)+(1/h^2)*rv1(j,1);
        end

        rv1(1,N)=rv1(1,N)*h^2/3;
        rv1(1,N-1)=rv1(1,N-1)+(1/h^2)*rv1(1,N);

        for i = N-1:-1:2
            rv1(1,i)=rv1(1,i)*h^2/4;
            rv1(1,i-1)=rv1(1,i-1)+(1/h^2)*rv1(1,i);
        end
        rv1(1,1)=rv1(1,1)*h^2/3;
        
        v(2:N,:)=v(2:N,:)+rv1;

    end
end