function A=FacA(N)

    %表示A
    A=sparse([1:N-1,N:(N-1)*(N-1),(N-1)*(N-1)+1:(N-1)*N,2:(N-1)*N,1:(N-1)*N-1,N:N*(N-1),1:(N-1)*(N-1)], ...
              [1:N-1,N:(N-1)*(N-1),(N-1)*(N-1)+1:(N-1)*N,1:(N-1)*N-1,2:(N-1)*N,1:(N-1)*(N-1),N:N*(N-1)], ...
        [linspace(3,3,N-1),linspace(4,4,(N-1)*(N-2)),linspace(3,3,N-1),linspace(-1,-1,N*(N-1)-1),linspace(-1,-1,N*(N-1)-1),linspace(-1,-1,(N-1)*(N-1)),linspace(-1,-1,(N-1)*(N-1))],N*(N-1),N*(N-1));
    A(N:N-1:(N-1)*(N-1)+1,N-1:N-1:(N-1)*(N-1))=zeros(N-1);
    A(N-1:N-1:(N-1)*(N-1),N:N-1:(N-1)*(N-1)+1)=zeros(N-1);

    if N<=128
        A=full(A);
    end

    % 对A做带状高斯消去分解
    for k = 1:(N-1)*(N-1)
        A(k+1:k+N-1,k) = A(k+1:k+N-1,k)./A(k,k);
		A(k+1:k+N-1,k+1:k+N-1) = A(k+1:k+N-1,k+1:k+N-1)-A(k+1:k+N-1,k)*A(k,k+1:k+N-1);
    end
    for k = (N-1)*(N-1)+1:N*(N-1)-1
        A(k+1:N*(N-1),k) = A(k+1:N*(N-1),k)./A(k,k);
        A(k+1:N*(N-1),k+1:N*(N-1)) = A(k+1:N*(N-1),k+1:N*(N-1))-A(k+1:N*(N-1),k)*A(k,k+1:N*(N-1));
    end
end