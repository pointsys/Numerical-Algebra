function b=Gauss(n,A,b)
	for k = 1:n-1
		A(k+1:n,k) = A(k+1:n,k)./A(k,k);
		A(k+1:n,k+1:n) = A(k+1:n,k+1:n)-A(k+1:n,k)*A(k,k+1:n);
	end
	%解下三角方程
	for k = 1:n-1
		b(k+1:n) = b(k+1:n)-b(k)*A(k+1:n,k)';
	end
	%解上三角方程
	for k = n:-1:2
		b(k) = b(k)/A(k,k);
		b(1:k-1) = b(1:k-1)-b(k)*A(1:k-1,k)';
	end
	b(1) = b(1)/A(1,1);
end