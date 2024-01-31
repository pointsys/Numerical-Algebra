A=diag(linspace(10,10,100))+diag(linspace(1,1,99),1)+diag(linspace(1,1,99),-1);
b=rand(1,100);
b=LDL(100,A,b);