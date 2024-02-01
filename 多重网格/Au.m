function Au=Au(m,u)
    Au=4*u;
    Au(2:m+1,2:m+1)=Au(2:m+1,2:m+1)-u(1:m,2:m+1);
    Au(2:m+1,2:m+1)=Au(2:m+1,2:m+1)-u(3:m+2,2:m+1);
    Au(2:m+1,2:m+1)=Au(2:m+1,2:m+1)-u(2:m+1,1:m);
    Au(2:m+1,2:m+1)=Au(2:m+1,2:m+1)-u(2:m+1,3:m+2);
end