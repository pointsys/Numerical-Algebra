for n=7:10

    eN=[];

    N=2^n;
    m=N-1;

    % 计算f
    x = linspace(0,1,N+1);
    y = linspace(0,1,N+1);
    [X,Y]=meshgrid(x,y);
    h=1/N;
    f = h^2*2*pi^2*sin(pi.*X).*sin(pi.*Y);
    f(m+2,:)=zeros(1,m+2);
    f(:,m+2)=zeros(m+2,1);

    %计算真解
    g = sin(pi*X).*sin(pi*Y);
    
    tic;
    u=Gauss_Seidel(zeros(m+2,m+2),5,f,m);
    r=f-Au(m,u);

    inte=0;

    R=1;
    e=[];

    while R>= 10^(-6)

        inte=inte+1;

        u=vcycle(u,r,N-1);
        r=f-Au(N-1,u);
        R=norm(r)/norm(f);
        e=[e,h*norm(u-g)];

    end
    toc;

    eN=[eN,e(inte)];

    plot(1:inte,e);
    title(strcat("N=",num2str(N)));
    
    saveas(gcf,strcat('N=',num2str(N),'.jpg'));

    clf;

end