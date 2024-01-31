n_rec=[];%用于存储迭代次数
for ep = [1,1e-1,1e-2,1e-3,1e-4,1e-5,1e-6,1e-7]
    %构建方程组
    A = [1,-1,0;-1,2,-1;0,-1,1]+ep*eye(3);
    b = [1,10,23]';
    a = mldivide(A,b);
    am = norm(a,2);
    x = [0,0,0]';

    n = 0;%用于记录迭代次数
    %进行改进Gauss-Seidel迭代
    while norm(x-a,2)/am > 1e-6
        %进行Gauss-Seidel迭代
        for k = 1:2
            r = b - A*x;
            e = r./[A(1,1),A(2,2),A(3,3)]';
            e(2) = e(2)-e(1)*A(2,1)/A(2,2);
            e(3) = e(3)-e(1)*A(3,1)/A(3,3)-e(2)*A(3,2)/A(3,3);
            x = x+e;
            n = n+1;
        end
        %对v_1方向进行修正
        r = b - A*x;
        alpha=sum(r)/sum(A,"all");
        x=x+alpha*[1,1,1]';
        n=n+1;
    end
    n_rec = [n_rec,n];
end