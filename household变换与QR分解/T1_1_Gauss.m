%表示A
A = diag(linspace(6,6,84))+diag(linspace(8,8,83),-1)+diag(linspace(1,1,83),1);

%表示b
b = linspace(15,15,84);
b(1) = 7;
b(84) = 14;

b=Gauss(84,A,b);

%求L2误差
b2=(b-1).^2;
L2_error=sum(b2);
L2_error = sqrt(L2_error);

%求L无穷误差
L_infin_error = max(abs(b-1));