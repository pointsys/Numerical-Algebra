function x=QR_LS_house(m,n,A,x)
%用Household方法做QR分解求LS问题的解
    %对矩阵A做QR分解
    Q=diag(linspace(1,1,m));
    for j = 1:n
        if j <m
            [v,b]=house(A(j:m,j)');
            Q1=diag(linspace(1,1,m));
            Q1(j:m,j:m)=Q1(j:m,j:m)-b.*v'*v;
            Q=Q*Q1;
            A(j:m,j:n)=Q1(j:m,j:m)*A(j:m,j:n);
        end
    end

    x=x*Q;
    
    %解上三角方程
    for k = n:-1:2
        x(k) = x(k)/A(k,k);
        x(1:k-1) = x(1:k-1)-x(k)*A(1:k-1,k)';
    end
    x(1) = x(1)/A(1,1);
    x=x(1:n);
end