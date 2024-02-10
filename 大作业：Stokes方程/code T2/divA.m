function [u,v]=divA(N,A,p,f,g)
    h=1/N;
    n=N*(N-1);
    
    f=f-Bup(N,p);
    f=f.*h^2;
    F=[];
    for i = 1:N
        F=[F,f(i,:)];
    end
    
    g=g-Bvp(N,p);
    g=g.*h^2;
    G=[];
    for j = 1:N
        G=[G,g(:,j)'];
    end

    % F=A\F';%*
    % G=A\G';%*
    % F=F';%*
    % G=G';%*
    
    %解下三角方程
	for k = 1:n-1
		F(k+1:n) = F(k+1:n)-F(k)*A(k+1:n,k)';
	end
	%解上三角方程
	for k = n:-1:2
		F(k) = F(k)/A(k,k);
		F(1:k-1) = F(1:k-1)-F(k)*A(1:k-1,k)';
	end
	F(1) = F(1)/A(1,1);

    %解下三角方程
	for k = 1:n-1
		G(k+1:n) = G(k+1:n)-G(k)*A(k+1:n,k)';
	end
	%解上三角方程
	for k = n:-1:2
		G(k) = G(k)/A(k,k);
		G(1:k-1) = G(1:k-1)-G(k)*A(1:k-1,k)';
	end
	G(1) = G(1)/A(1,1);

    u=[0,F(1:N-1),0];
    for i = 2:N
        u=[u;0,F((i-1)*(N-1)+1:i*(N-1)),0];
    end

    v=[0,G(1:N-1),0]';
    for i = 2:N
        v=[v,[0,G((i-1)*(N-1)+1:i*(N-1)),0]'];
    end
    
end