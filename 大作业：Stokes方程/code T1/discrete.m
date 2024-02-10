function [f,g,u_true,v_true,p_true]=discrete(N)
    h=1/N;
    x=linspace(0,1,N+1);
    x_u=x;
    y_u=x(1:N)+1/2*h;
    x_v=x(1:N)+1/2*h;
    y_v=x;
    x_p=x(1:N)+1/2*h;
    y_p=x_p;

    [X_u,Y_u]=meshgrid(x_u,y_u);
    [X_v,Y_v]=meshgrid(x_v,y_v);
    [X_p,Y_p]=meshgrid(x_p,y_p);

    u_true=(1-cos(2*pi*X_u)).*sin(2*pi*Y_u);
    v_true=-(1-cos(2*pi*Y_v)).*sin(2*pi*X_v);
    p_true=1/3*(X_p.^3)-1/12;

    f=-4*pi^2*(2*cos(2*pi*X_u(:,2:N))-1).*sin(2*pi*Y_u(:,2:N))+X_u(:,2:N).^2;
    f(1,:)=f(1,:)-(1/h)*2*pi*(1.-cos(2*pi*x_u(2:N)));
    f(N,:)=f(N,:)+(1/h)*2*pi*(1.-cos(2*pi*x_u(2:N)));
    g= 4*pi^2*(2*cos(2*pi*Y_v(2:N,:))-1).*sin(2*pi*X_v(2:N,:));
    g(:,1)=g(:,1)+1/h*(-2*pi*(cos(2*pi*y_v(2:N)') - 1));
    g(:,N)=g(:,N)+1/h*(2*pi*(cos(2*pi*y_v(2:N)') - 1));


end