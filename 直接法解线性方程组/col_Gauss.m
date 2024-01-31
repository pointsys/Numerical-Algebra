function b=col_Gauss(n,A,b)
	%对A进行列主元Gauss消去
P=diag(linspace(1,1,n)); %用于记录所作交换
for k = 1:n-1
    temp_max=abs(A(k,k)); %记录列元素的最大值
    temp_pos=k; %记录列元素最大值的位置
    for m = k+1:n
        if abs(A(m,k))>temp_max
            temp_max=abs(A(m,k));
            temp_pos=m;
        end
    end
	if temp_pos~=k %做交换
		l=A(k,:);
		A(k,:)=A(temp_pos,:);
		A(temp_pos,:)=l;
		l=P(k,:); %更新置换矩阵
		P(k,:)=P(temp_pos,:);
		P(temp_pos,:)=l;
	end
	A(k+1:n,k) = A(k+1:n,k)./A(k,k); %Gauss消去
	A(k+1:n,k+1:n) = A(k+1:n,k+1:n)-A(k+1:n,k)*A(k,k+1:n);
end

	%求Pb
	b=b*P';
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