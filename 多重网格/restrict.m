function F=restrict(m,f)
    F=f(3:2:m,2:2:m-1)+f(4:2:m+1,2:2:m-1)+f(4:2:m+1,3:2:m)+f(3:2:m,4:2:m+1)+f(2:2:m-1,4:2:m+1)+f(2:2:m-1,3:2:m);
    F=1/2*F;
    F=F+f(3:2:m,3:2:m);
    F=[zeros(1,(m+1)/2-1);F;zeros(1,(m+1)/2-1)];
    F=[zeros((m+1)/2+1,1),F,zeros((m+1)/2+1,1)];
end