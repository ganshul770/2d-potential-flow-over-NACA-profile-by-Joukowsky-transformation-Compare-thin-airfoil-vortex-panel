clear all;
close all;
clc;

%%
%joukowsky transformation
L=1;xmc=0.4*L;ymc=0.02*L;c=L/4;tm=0.12*c;ro=1.225;U=100;alpha=[0:pi/180:12*pi/180];n=160;
syms H
eqn=ymc==L*sqrt(0.25*(1+0.0625/H^2)-(xmc/L)^2)-0.125*L/H;
s1=L*double(solve(eqn,H));
n0=s1/2;
% n0=0;
g0=0.77*tm*c/L;
% g0=0;
delta=atan(n0/g0);
m=g0/cos(delta);
R=sqrt(m^2+c^2-2*m*c*cos(delta));
theta=2*pi/(n-1);
theta=0:theta:2*pi;
T0=g0+i*n0;
T=R*exp(i*theta)+T0;
Z=T+c^2./T;
x=(real(Z)+0.5*L); %0.5*L to shift airofoil to zero
y=(imag(Z));

for i=1:n-1
    dx(i)=x(i+1)-x(i);
    dy(i)=y(i+1)-y(i);
    x_cp(i)=(x(i+1)+x(i))/2;
    y_cp(i)=(y(i+1)+y(i))/2;
    l(i)=sqrt((x(i+1)-x(i))^2+(y(i+1)-y(i))^2);
end

%%
for o=1:length(alpha)
    [gamma,T_vel]=naca(n,x,y,l,x_cp,y_cp,alpha,o,U);

   
    %%
    %cl cm_le and cm_qc by panel method
    cl=0;
    cm_le=0;
    for j=1:n-1
        cl=cl+l(j)*(gamma(j)+gamma(j+1))/(U*L);
        cm_le=cm_le-l(j)*((2*x(j)*gamma(j)+x(j)*gamma(j+1)+x(j+1)*gamma(j)+2*x(j+1)*gamma(j+1))*cos(alpha(o))+(2*y(j)*gamma(j)+y(j)*gamma(j+1)+y(j+1)*gamma(j)+2*y(j+1)*gamma(j+1))*sin(alpha(o)))/(3*U*L^2);
    end
    cm_qc=cm_le+cl*(0.25);


    %%
    %cl_c cm_le_c and cm_qc_c by complex potential and by thin airofoil theory

    cl_c=2*pi*(sin(alpha(o))+n0*cos(alpha(o))/sqrt(R^2-n0^2))/(1+g0/(sqrt(R^2-n0^2)-g0));
    cm_le_c=-(0.25*pi*sin(2*alpha(o))*(R^2-n0^2-g0^2).^2/(R^2-n0^2).^2)-(0.25*cl_c*((0.25-g0)*cos(alpha(o))-n0*sin(alpha(o)))*(sqrt(R^2-n0^2)-g0)/(R^2-n0^2))-(0.0625*cl_c*cos(alpha(o))*(sqrt(R^2-n0^2)-g0)/(R^2-n0^2));
    cm_qc_c=-0.25*pi*sin(2*alpha(o))*(R^2-n0^2-g0^2).^2/(R^2-n0^2).^2-0.25*cl_c*((0.25-g0)*cos(alpha(o))-n0*sin(alpha(o)))*(sqrt(R^2-n0^2)-g0)/(R^2-n0^2);

    %%
    %camber
    Y=L*sqrt(0.25*(1+0.0625/s1^2)-((x(1:1+0.5*n)-0.5*L)./L).^2)-0.125*L/s1;
    %cl_t cm_le_t and cm_qc_t by thin airofoil theory
    for i=1:80
        dyc_dx(i)=-(Y(i+1)-Y(i))/(x(i+1)-x(i));
    end
%     cl_t=2*pi*(alpha(o)-(1/pi)*trapz(theta(1:80),dyc_dx.*(1-cos(theta(1:80)))));
    cl_t=2*pi*(sqrt(R^2-n0^2)-g0)*(alpha(o)-atan(-n0/sqrt(R^2-n0^2)))/sqrt(R^2-n0^2);
    cm_le_t=-0.25*cl_t+0.5*trapz(theta(1:80),dyc_dx.*(cos(2*theta(1:80))-cos(theta(1:80))));
    cm_qc_t=cm_le_t+0.25*cl_t/L;
    figure(1)
    plot(180*alpha(o)/pi,cl,'-o');
    hold on;
    plot(180*alpha(o)/pi,cl_c,'-p');
    hold on;
    plot(180*alpha(o)/pi,cl_t,'-s');
    hold on;
    xlabel('---- alpha(in degree) ---->');
    ylabel('---- lift coefficient ---->');
    legend('panel','joukowsky','thin airfoil');

    figure(2)
    plot(180*alpha(o)/pi,cm_le,'-o');
    hold on;
    plot(180*alpha(o)/pi,cm_le_c,'-p');
    hold on;
    plot(180*alpha(o)/pi,cm_le_t,'-s');
    hold on;
    xlabel('---- alpha(in degree) ---->');
    ylabel('---- leading edge moment coefficient ---->');
    legend('panel','joukowsky','thin airfoil');
    axis([0 12 -3 3]);
    
    figure(3)
    plot(180*alpha(o)/pi,cm_qc,'-o');
    hold on;
    plot(180*alpha(o)/pi,cm_qc_c,'-p');
    hold on;
    plot(180*alpha(o)/pi,cm_qc_t,'-s');
    hold on;
    xlabel('---- alpha(in degree) ---->');
    ylabel('---- quater chord moment coefficient ---->');
    legend('panel','joukowsky','thin airfoil');
    axis([0 12 -3 3]);
end
%%
% joukowsky profile
% alpha=pi*input('enter the value of alpha = ')/180;
alpha=0*pi/180;o=1;
[gamma,T_vel]=naca(n,x,y,l,x_cp,y_cp,alpha,o,U);
%camber
Y=L*sqrt(0.25*(1+0.0625/s1^2)-((x(1:1+0.5*n)-0.5*L)./L).^2)-0.125*L/s1;
figure(4)
plot(x,y);
hold on;
plot(x(1:1+0.5*n),Y);
xlabel('---- x/c ---->');
ylabel('---- y/c ---->');
title('profile with camber');

%%
    %cl cm_le and cm_qc by panel method
    cl=0;
    cm_le=0;
    for j=1:n-1
        cl=cl+l(j)*(gamma(j)+gamma(j+1))/(U*L);
        cm_le=cm_le-l(j)*((2*x(j)*gamma(j)+x(j)*gamma(j+1)+x(j+1)*gamma(j)+2*x(j+1)*gamma(j+1))*cos(alpha(o))+(2*y(j)*gamma(j)+y(j)*gamma(j+1)+y(j+1)*gamma(j)+2*y(j+1)*gamma(j+1))*sin(alpha(o)))/(3*U*L^2);
    end
    cm_qc=cm_le+cl*(0.25);


    %%
    %cl_c cm_le_c and cm_qc_c by complex potential and by thin airofoil theory

    cl_c=2*pi*(sin(alpha(o))+n0*cos(alpha(o))/sqrt(R^2-n0^2))/(1+g0/(sqrt(R^2-n0^2)-g0));
    cm_le_c=-(0.25*pi*sin(2*alpha(o))*(R^2-n0^2-g0^2).^2/(R^2-n0^2).^2)-(0.25*cl_c*((0.25-g0)*cos(alpha(o))-n0*sin(alpha(o)))*(sqrt(R^2-n0^2)-g0)/(R^2-n0^2))-(0.0625*cl_c*cos(alpha(o))*(sqrt(R^2-n0^2)-g0)/(R^2-n0^2));
    cm_qc_c=-0.25*pi*sin(2*alpha(o))*(R^2-n0^2-g0^2).^2/(R^2-n0^2).^2-0.25*cl_c*((0.25-g0)*cos(alpha(o))-n0*sin(alpha(o)))*(sqrt(R^2-n0^2)-g0)/(R^2-n0^2);

    %%
%     cl_t cm_le_t and cm_qc_t by thin airofoil theory
    for i=1:80
        dyc_dx(i)=-(Y(i+1)-Y(i))/(x(i+1)-x(i));
    end
%     cl_t=2*pi*(alpha(o)-(1/pi)*trapz(theta(1:80),dyc_dx.*(1-cos(theta(1:80)))));
    cm_le_t=-0.25*cl_t+0.5*trapz(theta(1:80),dyc_dx.*(cos(2*theta(1:80))-cos(theta(1:80))));
    cm_qc_t=cm_le_t+0.25*cl_t/L;
    cl_t=2*pi*(sqrt(R^2-n0^2)-g0)*(alpha(o)-atan(-n0/sqrt(R^2-n0^2)))/sqrt(R^2-n0^2);

%%
%cp variation
cp=1-(T_vel.^2./U^2);
figure(5)
plot(x_cp,cp(1:n-1));
xlabel('---- x/c ---->');
ylabel('---- cp ---->');
title("cp on upper and lower surface");
set(gca, 'YDir','reverse');

%%
% % cp variation by complex potential method
% for j=1:n-1
%     T_cp(j)=(T(j+1)+T(j))/2;
% end
% w_z=(exp(-i*alpha)+(2*i*(sqrt(R^2-n0^2)*sin(alpha)+n0*cos(alpha))./(T_cp-T0))-(R^2*exp(i*alpha)./(T_cp-T0).^2))./(1-c^2./T_cp.^2);
% cp_c=(1-w_z.^2);
% figure(4)
% plot(x_cp,abs(cp_c(1:n-1)));
% set(gca, 'YDir','reverse');

%%
a=[x(0.5*n-1),0,x(0.5*n+1)];
b=[y(0.5*n-1),0,y(0.5*n+1)];
vt=zeros(3,1);
for i=1:3
    h=zeros(2,1);
    [h]=global_vel(x,y,l,a(i)+10^(-4),b(i)+10^(-4),h,gamma,n,alpha,U);

    vt(i)=[(x(0.5*n-1+i)-a(i))/sqrt((x(0.5*n-1+i)-a(i))^2+(y(0.5*n-1+i)-b(i))^2),(y(0.5*n-1+i)-b(i))/sqrt((x(0.5*n-1+i)-a(i))^2+(y(0.5*n-1+i)-b(i))^2)]*h;
end
d_u_vt=-(vt(3)-vt(2))/(a(3)-a(2));
d_l_vt=(vt(2)-vt(1))/(a(2)-a(1));


%%
%search stagnation point
vt_le=vt(2);
N_R=0;
u_vt_n=0;
l_vt_n=0;
while (u_vt_n<10^(-3)*U || l_vt_n<10^(-3)*U)
    if (vt_le<0)%search upper
        N_R=N_R-vt_le/d_u_vt;
    else%search lower
        N_R=N_R-vt_le/d_l_vt;
    end
    xl_g=0;yl_g=0;xu_g=0;yu_g=0;
    %     [xl_n,yl_n,xu_n,yu_n]=surface_coordinate(tm,xmc_c,ymc_c,L,N_R);
    geta=acos(N_R/(R+c^2/R));
    xl_n=N_R;
    xu_n=N_R;
    yu_n=(R-c^2/R)*sin(geta);
    yl_n=(R-c^2/R)*sin(pi+geta);

    h_u=zeros(2,1);
    h_l=zeros(2,1);
    [h_u]=global_vel(x,y,l,xu_n,yu_n,h_u,gamma,n,alpha,U);
    [h_l]=global_vel(x,y,l,xl_n,yl_n,h_l,gamma,n,alpha,U);
    xl_n=[xl_g,xl_n];yl_n=[yl_g,yl_n];xu_n=[xu_g,xu_n];yu_n=[yu_g,yu_n];
    u_vt_n=[(xu_n(2)-xu_n(1))/sqrt((xu_n(2)-xu_n(1))^2+(yu_n(2)-yu_n(1))^2),(yu_n(2)-yu_n(1))/sqrt((xu_n(2)-xu_n(1))^2+(yu_n(2)-yu_n(1))^2)]*h_u;
    l_vt_n=[(xl_n(2)-xl_n(1))/sqrt((xl_n(2)-xl_n(1))^2+(yl_n(2)-yl_n(1))^2),(yl_n(2)-yl_n(1))/sqrt((xl_n(2)-xl_n(1))^2+(yl_n(2)-yl_n(1))^2)]*h_l;
    xl_g=xl_n(2);yl_g=yl_n(2);xu_g=xu_n(2);yu_g=yu_n(2);
end

if u_vt_n<l_vt_n
    xs=xu_g;
    ys=yu_g;
    v_s=h_u;
else
    xs=xl_g;
    ys=yl_g;
    v_s=h_l;
end

%%
%plot
%from trailing edge
[t_T,s_T]=sl_by_point(x,y,l,gamma,n,alpha,U,1.001,0,0.001,25);
%from leading edge
[t_L,s_L]=sl_by_point(x,y,l,gamma,n,alpha,U,xs,ys+0.006,-0.001,25);

figure(6)
plot(x,y);
xlabel('---- x/c ---->');
ylabel('---- y/c ---->');
title("streamline");
hold on;
plot(t_T,s_T);
hold on;
plot(t_L,s_L);
hold on;
%other streamlines
for i=-21:2:21
    [m,o]=sl_by_point(x,y,l,gamma,n,alpha,U,t_L(end),s_L(end)+0.03*i,0.001,60);
    plot(m,o);
    hold on;
end
