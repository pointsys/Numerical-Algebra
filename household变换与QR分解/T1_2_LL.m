A=diag(linspace(10,10,100))+diag(linspace(1,1,99),1)+diag(linspace(1,1,99),-1);
b = linspace(12,12,100);
b(1) = 11;
b(100) = 11;
x=LL(100,A,b);

%求L2误差
b2=(x-1).^2;
L2_error=sum(b2);
L2_error = sqrt(L2_error);

%求L无穷误差
L_infin_error = max(abs(x-1));