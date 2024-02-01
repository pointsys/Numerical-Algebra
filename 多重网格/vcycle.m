function e=vcycle(e0,r,m)
    if m==1
        e=e0+r/4;
    else

        r1=restrict(m,r);%限制残量
        e1=Gauss_Seidel(zeros((m+1)/2+1,(m+1)/2+1),5,r1,(m+1)/2-1);%以0为初值，对残量方程进行迭代
        r2=r1-Au((m+1)/2-1,e1);%新的残量
    
        e2=vcycle(e1,r2,(m+1)/2-1);%粗网格上误差，即细网格残量方程的近似解
        
        e3=lift(e2,(m+1)/2-1);%细网格上误差进行提升

        e4=Gauss_Seidel(e3,5,r,m);%对残量方程再磨光

        e=e0+e4;%更新误差

    end

end