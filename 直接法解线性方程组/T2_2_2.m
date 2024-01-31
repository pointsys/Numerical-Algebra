for i = 1:40
	for j = 1:40
		A(i,j)=1/(i+j-1);
	end
	b(i)=sum(A(i,:));
end
b=LDL(40,A,b);