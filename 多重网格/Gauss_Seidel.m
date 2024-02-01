function u=Gauss_Seidel(u,v1,f,m)
%进行v1步Gauss-Seidel迭代

    for k=1:v1

        r=f-Au(m,u);

        for i=2:m
            for j = 2:m+1
                r(i,j)=r(i,j)/4;
                r(i+1,j)=r(i+1,j)+r(i,j);
                if j ~=m+1
                    r(i,j+1)=r(i,j+1)+r(i,j);
                end
            end
        end

        for j = 2:m
            r(m+1,j)=r(m+1,j)/4;
            r(m+1,j+1)=r(m+1,j+1)+r(m+1,j);
        end
        r(m+1,m+1)=r(m+1,m+1)/4;
        
        u=u+r;
    end
    
end