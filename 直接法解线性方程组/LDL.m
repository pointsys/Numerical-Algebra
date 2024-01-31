function b = LDL(n,A,b)
	for j = 1:n
		for i = 1:j-1
			v(i)=A(j,i)*A(i,i);
		end
		if j == 1
			A(j+1:n,j)=A(j+1:n,j)/A(j,j);
		elseif j == n
			A(n,n)=A(n,n)-A(n,1:n-1)*v(1:n-1)';
		else
			A(j,j)=A(j,j)-A(j,1:j-1)*v(1:j-1)';
			A(j+1:n,j)=(A(j+1:n,j)-A(j+1:n,1:j-1)*v(1:j-1)')/A(j,j);
		end
	end
	%解下三角方程
	for k = 1:n-1
		b(k+1:n) = b(k+1:n)-b(k)*A(k+1:n,k)';
	end
	%解对角方程
	for i =1:n
		b(i)=b(i)/A(i,i);
	end
	%解上三角方程
	for k = n:-1:2
		b(1:k-1) = b(1:k-1)-b(k)*A(k,1:k-1);
	end
end	