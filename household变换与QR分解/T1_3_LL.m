for i = 1:40
	for j = 1:40
		A(i,j)=1/(i+j-1);
	end
	b(i)=sum(A(i,:));
end
x=LL(40,A,b);

%求L2误差
x2=(x-1).^2;
L2_error=sum(x2);
L2_error = sqrt(L2_error);

%求L无穷误差
L_infin_error = max(abs(x-1));