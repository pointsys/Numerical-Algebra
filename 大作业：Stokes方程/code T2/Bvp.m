function Bvp=Bvp(N,p)
    h=1/N;
    Bvp=1/h*(p(2:N,1:N)-p(1:N-1,1:N));
end