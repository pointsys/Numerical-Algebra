function Bup=Bup(N,p)
    h=1/N;
    Bup=1/h*(p(1:N,2:N)-p(1:N,1:N-1));
end