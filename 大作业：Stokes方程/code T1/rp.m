function rp=rp(N,u,v,d)
    h=1/N;
    rp=d+(u(1:N,2:N+1)-u(1:N,1:N)+v(2:N+1,1:N)-v(1:N,1:N))/h;
end