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
% 对A做带状高斯消去分解
    for i = 1:N-2
        for k = 1:N-1
            m=(i-1)*(N-1)+k;
            t=(i+1)*(N-1);
            A(m+1:t,m)=A(m+1:t,m)/A(m,m);
            A(m+1:t,m+1:t)=A(m+1:t,m+1:t)-A(m+1:t,m)*A(m,m+1:t);
        end
    end
    for k = (N-1)*(N-2)+1:n-1
        A(k+1:n,k)=A(k+1:n,k)/A(k,k);
        A(k+1:n,k+1:n)=A(k+1:n,k+1:n)-A(k+1:n,k)*A(k,k+1:n);
    end

% 计算得解
    for k = 2:n
        f(k:n)=f(k:n)-A(k:n,k-1)'*f(k-1);
    end
    for k = n-1:-1:1
        f(k+1)=f(k+1)/A(k+1,k+1);
        f(1:k)=f(1:k)-A(1:k,k+1)'*f(k+1);
    end
    f(1)=f(1)/A(1,1);
    toc;
end