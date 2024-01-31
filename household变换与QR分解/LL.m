function b=LL(n,A,b)
	for i = 1:n-1
		A(i,i)=sqrt(A(i,i));
		A(i+1:n,i)=A(i+1:n,i)/A(i,i);
		for j =i+1:n
			A(j:n,j)=A(j:n,j)-A(j:n,i)*A(j,i);
		end
	end		
	A(n,n)=sqrt(A(n,n));
	%解下三角方程
	for k = 1:n-1
		b(k) = b(k)/A(k,k);
		b(k+1:n) = b(k+1:n)-b(k)*A(k+1:n,k)';
	end
	b(n) = b(n)/A(n,n);
	%解上三角方程
	for k = n:-1:2
		b(k) = b(k)/A(k,k);
		b(1:k-1) = b(1:k-1)-b(k)*A(k,1:k-1);
	end
	b(1) = b(1)/A(1,1);
end