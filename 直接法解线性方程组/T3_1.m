for N = [32,64,128]
    n = (N-1)*(N-1);
    
%表示A
    A = eye((N-1)*(N-1))*4+diag(linspace(-1,-1,(N-1)*(N-1)-1),1)+diag(linspace(-1,-1,(N-1)*(N-1)-1),-1);
    for i = 1:N-2
        A(i*(N-1)+1,i*(N-1))=0;
        A(i*(N-1),i*(N-1)+1)=0;
    end 
    for i = 1:(N-1)*(N-2)
        A(N+i-1,i)=-1;
        A(i,N+i-1)=-1;
    end

% 计算f
    x = linspace(0,1,N+1);
    x = x(2:N);
    y = linspace(0,1,N+1);
    y = y(2:N);
    f = 0;
    for t = y
        f = [f,sin(pi*x)*sin(pi*t)];
    end
    f = f(2:(N-1)*(N-1)+1);
        
    tic;
%LDL'分解
    v=0;
    for j = 1:n
        for i = 1:j-1
            v(i)=A(j,i)*A(i,i);
        end
        A(j,j)=A(j,j)-A(j,1:j-1)*v(1:j-1)';
        A(j+1:n,j)=(A(j+1:n,j)-A(j+1:n,1:j-1)*v(1:j-1)')/A(j,j);
    end
%解方程
    for k = 2:n 
        f(k:n)=f(k:n)-A(k:n,k-1)'*f(k-1);
    end
    for k = 1:n
        f(k)=f(k)/A(k,k);
    end
    for k = n-1:-1:1
        f(1:k)=f(1:k)-A(k+1,1:k)*f(k+1);
    end
    toc; 
end