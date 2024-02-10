function Av=Av(N,v,p)
    h=1/N;
    v=[zeros(N+1,1),v,zeros(N+1,1)];%添加ghost point
    v(:,1)=v(:,2);
    v(:,N+2)=v(:,N+1);
    Av=(1/h^2)*(4*v(2:N,2:N+1)-v(1:N-1,2:N+1)-v(3:N+1,2:N+1)-v(2:N,1:N)-v(2:N,3:N+2));
    Av=Av+(p(2:N,1:N)-p(1:N-1,1:N))/h;
end