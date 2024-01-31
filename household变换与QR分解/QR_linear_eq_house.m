function x=QR_linear_eq_house(n,A,x)
%用Household方法做QR分解求解线性方程组
    %对矩阵A做QR分解
    Q=diag(linspace(1,1,n));
    for j = 1:n-1
        [v,b]=house(A(j:n,j)');
        Q1=diag(linspace(1,1,n));
        Q1(j:n,j:n)=Q1(j:n,j:n)-b.*v'*v;
        Q=Q*Q1;
        A(j:n,j:n)=Q1(j:n,j:n)*A(j:n,j:n);
    end

    x=x*Q;

    %解上三角方程
    for k = n:-1:2
        x(k) = x(k)/A(k,k);
        x(1:k-1) = x(1:k-1)-x(k)*A(1:k-1,k)';
    end
    x(1) = x(1)/A(1,1);
end