function [gamma,T_vel]=naca(n,x,y,l,x_cp,y_cp,alpha,o,U)
    A=zeros(n,n);
    Tn=zeros(n,n);
    B=zeros(n,1);
    v_in_t=zeros(n,1);
    for i=1:n-1
        for j=1:n-1
            [P]=influence_matrix_vortex(x(j),x(j+1),y(j),y(j+1),l(j),x_cp(i),y_cp(i));
            A(i,j)=A(i,j)+(-P(1,1)*(y(i+1)-y(i))/l(i)+P(2,1)*(x(i+1)-x(i))/l(i));
            A(i,j+1)=A(i,j+1)+(-P(1,2)*(y(i+1)-y(i))/l(i)+P(2,2)*(x(i+1)-x(i))/l(i));
        end
    end
    for i=1:n-1
        for j=1:n-1
            [P]=influence_matrix_vortex(x(j),x(j+1),y(j),y(j+1),l(j),x_cp(i)+0.001*(-(y(i+1)-y(i))/l(i)),y_cp(i)+0.001*((x(i+1)-x(i))/l(i)));
            Tn(i,j)=Tn(i,j)+(P(1,1)*(x(i+1)-x(i))/l(i)+P(2,1)*(y(i+1)-y(i))/l(i));
            Tn(i,j+1)=Tn(i,j+1)+(P(1,2)*(x(i+1)-x(i))/l(i)+P(2,2)*(y(i+1)-y(i))/l(i));
        end
    end

    A(n,[1,end])=1;
    Tn(n,[1,end])=1;

    for i=1:n-1
        B(i,1)=U*((y(i+1)-y(i)).*cos(alpha(o))-(x(i+1)-x(i)).*sin(alpha(o)))./l(i);
        v_in_t(i,1)=U*((x(i+1)-x(i))*cos(alpha(o))+(y(i+1)-y(i))*sin(alpha(o)))/l(i);
    end
    gamma=A\B;
    T_vel=Tn*gamma+v_in_t;