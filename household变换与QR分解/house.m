function [v,b]=house(x)
%对给定向量做Household变换，使得该向量除了首个分量外都变为0
    n=length(x);
    eta=max(abs(x));
    x=x/eta;
    sigma=x(2:n)*x(2:n)';
    v=0;
    v(2:n)=x(2:n);
    if sigma==0
        b=0;
    else
        a=sqrt(x(1)^2+sigma);
        if x(1)<=0
            v(1)=x(1)-a;
        else
            v(1)=-sigma/(x(1)+a);
        end
        b=2*v(1)^2/(sigma+v(1)^2);
        v=v/v(1);
    end
end
