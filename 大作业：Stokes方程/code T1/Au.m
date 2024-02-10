function Au=Au(N,u,p)
    h=1/N;
    u=[zeros(1,N+1);u;zeros(1,N+1)];%添加ghost point
    u(1,:)=u(2,:);
    u(N+2,:)=u(N+1,:);
    Au=4*u(2:N+1,2:N)-u(2:N+1,1:N-1)-u(2:N+1,3:N+1)-u(1:N,2:N)-u(3:N+2,2:N);
    Au=1/h^2*Au+1/h*(p(1:N,2:N)-p(1:N,1:N-1));
end