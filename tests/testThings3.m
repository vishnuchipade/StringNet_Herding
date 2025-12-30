%%
% ------------------------------
% Author: Vishnu S. Chipade
% Email: vishnuc@umich.edu
% ------------------------------

%test k-means clustering with size constraints on the clusters depending on
%the cardinality of the cluster

circCos=cos(0:pi/100:2*pi);
circSin=sin(0:pi/100:2*pi);
markSize=4;

colors={[1,0,0],[0,1,0],[0,0,1],[0,0,0],[0.5,0.5,0.5],[0,0.5,1],[1,0.5,1],[0.5,0.5,1],[1,.1,0.5],[0.1,0.9,0.4],[0.9,0.4,0.5],[0,0.5,0.4]};


NA=10;
if (0)
    rA=rand(2,NA);
    vA=rand(2,NA);
end
XA=[rA]%;vA];
XA=XA';
NC=2;
[idx,C,sumd,D] =kmeans(XA,NC);

%plot the points
figure
hold all;
plot(XA(:,1),XA(:,2),'b*')
%plot the centers and circles around the center with larges

for k=1:NC
    %center of the cluster
    plot(C(k,1),C(k,2),'o','color',colors{k},'markersize',markSize)
    ind=find(idx==k);
    for ki=1:length(ind)
        Dk(ki)=norm(XA(ind(ki),1:2)-C(k,1:2));
    end
    Rk=max(Dk);
    plot(XA(ind,1),XA(ind,2),'o','color',colors{k})
    plot(C(k,1)+Rk*circCos,C(k,2)+Rk*circSin,'--','color',colors{k});
    text(C(k,1),C(k,2),['$C_',num2str(k),'$'])
end


figure
gscatter(XA(:,1),XA(:,2),idx)

%using dbscan
if(0)
    idx2 = dbscan(XA,1,1)
    
    figure
    gscatter(XA(:,1),XA(:,2),idx2)
end



%test distance travelled at extreme speeds/accelerations for modified
%double integrator
%%
u_ma=5;
Cd=0.1;
v_ma=sqrt(u_ma/Cd);
s_ma=log(2)/(2*Cd)
s(1)=0;
t(1)=0;
dt=0.001;
v(1)=v_ma;

i=1;
while v(i)>0
    s(i+1)=s(i)+dt*v(i);
    v(i+1)=v(i)+dt*(-u_ma-Cd*v(i)^2);
    i=i+1;
end
s_f=s(i)
v_f=v(i)

%%
%Double integrator along circular arc
dt=0.001;
v(1)=0;
C_d=0.1;
u=10;
rho=1;
v_m=sqrt(0.5*(sqrt(C_d^2*rho^4+4*rho^2*u^2)-C_d*rho^2))
for i=1:10000
    v(i+1)=v(i)+dt*(-C_d*v(i)^2+(u-v(i)^4/(u*rho^2)));
end
figure
plot(v)

%%
%plot x-x^alpha for x=[0,x_m]
x_m=10;
alpha=0.5;
C_d=0.01;
kr=1.7;
x=0:x_m/100:x_m;
figure
plot(x,C_d*x-kr*x.^(alpha))


%%
%Check Lyapunov function and its derivative for Double integrator with bounded control


SetPlotDefaults
global C_d dt
C_d=.1;
dt=0.001;
NSim=200000;
u_max=2;
%v_max=u_max/C_d;
v_max=sqrt(u_max/C_d);
eta_max=v_max/3;
u_max2=2*C_d*eta_max^2;
u_max1=u_max-u_max2;
%v_max=sqrt(u_max/C_d);

r0=10/sqrt(2)*[1,1]';
xi=[0,0]';
theta0=atan2(r0(2)-xi(2),r0(1)-xi(1))+pi/2*rand(1);
theta0=pi;
v0=eta_max+(v_max-eta_max)*rand(1);
v0=v_max;
X(:,1)=[r0; v0*cos(theta0);v0*sin(theta0)];

kr=.2;
alpha=6.01;
kr_opt=(-alpha+sqrt(alpha^2+4*(u_max^3/C_d^2)))/(2*u_max/C_d^2);

%kr=kr_opt;

A=zeros(4);
A(1:2,3:4)=eye(2);
A(3:4,1:2)=-kr*eye(2);
A(3:4,3:4)=-C_d*eye(2);
P=lyap(A',eye(4));
lambda_min=min(eig(P));


%Symbolic expressions for Lyapunov function V=x'Px
syms x1 x2 x3 x4 p1 p2 p12 sigma_r1 sigma_r2 cd k
As=[zeros(2),eye(2);-k*eye(2),-cd*eye(2)];
Ps=[p1*eye(2),p12*eye(2);p12*eye(2),p2*eye(2)];
Xs=[x1;x2;x3;x4];
Xs_dot=[x3;x4;sigma_r1-cd*x3;sigma_r2-cd*x4];
Vs=Xs'*Ps*Xs
Vs_dot=Xs_dot'*Ps*Xs+Xs'*Ps*Xs_dot


kv=1;
alpha_v=1;
alpha_r=alpha_v/(2-alpha_v);

R_max=u_max/kr;

% R_max1=((abs(C_d*v_max-kv*v_max^(alpha_v))+u_max)/kr)^(1/alpha_r)
% R_max2=u_max/C_d^2
% R_max3=(2*u_max/kr)^(1/alpha_r);

%Invariant Level set
p1=P(1,1);
p2=P(4,4);
p12=P(1,3);
Q0=1;
c=p1*(R_max)^2%+p2*v_max^2+2*p12*R_max*v_max
R_max2=(2*p2*u_max*v_max+Q0*v_max^2)/(2*p1*v_max-4*p12*u_max)
c0=p1*(R_max2)^2;%+p2*v_max^2+2*p12*R_max2*v_max

%R_maxV=sqrt(c/p1)

x=0:v_max/100:v_max;
cm1=0;
cm2=1050000;
y_max0=cm2;
y_max=cm2;
%Plot max V_dot
for j=1:1000
    cm=0+1000*j;
    vm=min(v_max,0.99999*sqrt(4*p1*cm/abs(4*p1*p2-4*p12^2)));
    x=0:vm/100:vm;
    for i=1:length(x)
        r_norm=(-2*p12*x(i)+sqrt((2*p12*x(i))^2-4*p1*(p2*x(i)^2-cm)))/(2*p1);
        y(i)=((2*p1-2*p12*C_d)*x(i)-2*p12*u_max)*r_norm-2*p2*u_max*x(i)-Q0*x(i)^2;
    end
    y_max(j)=max(y);
end
figure
plot(y_max)

while (abs(y_max)>1e-5)
    cm=(cm1+cm2)/2;
    vm=min(v_max,0.99999*sqrt(4*p1*cm/abs(4*p1*p2-4*p12^2)));
    vm=min(v_max,sqrt(cm/p2));
    x=0:vm/100:vm;
    for i=1:length(x)
        r_norm=(-2*p12*x(i)+sqrt((2*p12*x(i))^2-4*p1*(p2*x(i)^2-cm)))/(2*p1);
        y(i)=((2*p1-2*p12*C_d)*x(i)-2*p12*u_max)*r_norm-2*p2*u_max*x(i)-Q0*x(i)^2;
    end
    y_max=max(y);
    if y_max>0
        cm2=(cm1+cm2)/2;
    else
        cm1=(cm1+cm2)/2;
    end
    err=y_max-y_max0;
    y_max0=y_max;
end
vm=min(v_max,0.99999*sqrt(4*p1*cm/abs(4*p1*p2-4*p12^2)));
vm=min(v_max,0.99999*sqrt(cm/p2))
x=0:vm/100:vm;
for i=1:length(x)
    r_norm_cm(i)=(-2*p12*x(i)+sqrt((2*p12*x(i))^2-4*p1*(p2*x(i)^2-cm)))/(2*p1);
end
R_max_inv1=min(r_norm_cm)%(-2*p12*v_max+sqrt((2*p12*v_max)^2-4*p1*(p2*v_max^2-cm)))/(2*p1)
R_max_inv2=sqrt(cm/p1);

R_minV=(v_max^2/kr)^(1/3);
R_minV2=sqrt(v_max^2/(u_max))

v_star=(4*p12^2+c)/(2*Q0+4*p12*(2*p1-2*p12*C_d));
f=@(x) (2*p1-2*p12*C_d)*sqrt(4*p1*x)/(2*p1)-2*p2*u_max-sqrt((4*Q0*p12*u_max*sqrt(4*p1*x))/(p1))
cm_approx=fsolve(f,100)
c_unsat=p1*v_max^2/u_max+p2*v_max^2%+2*p12*v_max^2/sqrt(u_max)
cm
cm-c_unsat
%R_maxV=sqrt(cm/p1);


R_max_FT2=1/kr*abs(abs(C_d*v_max-kv*v_max^(alpha_v))+u_max);
R_max_FT1=1/kr*abs(abs(C_d*v_max-kv*v_max^(alpha_v))-u_max);

vel_norm(1)=norm(X(3:4,1));
t(1)=0;

for i=1:NSim
    if i>37000
        
    end
    xi=[eta_max,0]'*t(i);
    eta=[eta_max,0]';
    eta_dot=[0,0]';
    dr=X(1:2,i)-xi;
    dv=X(3:4,i)-eta;
    norm_dr=norm(dr);
    norm_v= norm(X(3:4,i));
    norm_dv= norm(dv);
    %linear asymptotic controller
    %U=-kr*(X(1:2,i)-xi);
    
    %finite-time controller
    if norm_dr>0
        h1=-kr*(X(1:2,i)-xi)*norm_dr^(alpha_r-1);
    else
        h1=[0,0]';
    end
    if norm_v>0
       % h2=-kv*(dv)*norm_dv^(alpha_v-1)+C_d*X(3:4,i);
        h2=-kv*(dv)*norm_dv^(alpha_v-1)+C_d*norm_v*X(3:4,i);
    else
        h2=[0,0]';
    end
    flag=0;
    if (flag)
        U=h1+h2;
        norm_U(i)=norm(U);
        U=min(u_max,norm_U(i))*U/norm_U(i);  %Saturation;
        
        flag_sat(i)=0;
        if norm_U(i)>u_max
            flag_sat(i)=1;
        end
    else
        norm_h1=norm(h1);
        norm_h2=norm(h2);
        U=eta_dot+min(u_max1,norm_h1)*h1/norm_h1+min(u_max2,norm_h2)*h2/norm_h2;
        norm_U(i)=norm(U);
        if norm_h1>u_max1 && norm_h2>u_max2
            flag_sat(i)=3;
        elseif norm_h1>u_max1
           flag_sat(i)=1; 
        elseif norm_h2>u_max2
            flag_sat(i)=2;
        else
            flag_sat(i)=0;
        end
    end
    
    %Dynamics with linear drag term
%     dv_dot=U-C_d*X(3:4,i)-eta_dot;
%     X(:,i+1)=modifiedDIDynamicsLin(X(:,i),U);
    
    %dynamics with quadratic term
    dv_dot=U-C_d*norm_v*X(3:4,i)-eta_dot;
    X(:,i+1)=modifiedDIDynamics(X(:,i),U);
    %Lyapunov function 1
    
    dist_norm(i)=norm_dr;
    dvr=(X(3:4,i)'*(dr));
    norm_v=norm(X(3:4,i));
    V1(i)=0.5*norm_dr^2;
    V1_dot(i)=X(3:4,i)'*(X(1:2,i)-xi);
    V1_ddot(i)=dv_dot'*dr+norm_dv^2;
    
    %Lyapunov function 2
    V2(i)=1/(alpha_r+1)*norm_dr^(alpha_r+1)+norm_dv^2/2;
    V2_dot(i)=dvr*norm_dr^(alpha_r-1);
    if alpha_r<1
        V2_ddot(i)=(dv_dot'*dr+norm_dv^2)*norm_dr^(alpha_r+1)+dvr^2*norm_dr^(alpha_r-3)/(alpha_r-1);
    else
        V2_ddot(i)= dv_dot'*dr+norm_dv^2;
    end
    
    %Lyapunov function 3
    V3(i)=1/(alpha_r+1)*norm_dr^(alpha_r+1)+0.5*norm_dv^2;
    V3_dot(i)=dvr*norm_dr^(alpha_r-1)+dv_dot'*dv;
    %     V3(i) = 6.25*norm_dr^2+105*norm_v^2+20*norm_dr*norm_v;
    %     V3_dot(i)=-norm_dr^2-norm_v^2;
    dX=X(:,i)-[xi;eta];
    %     V3(i)=dX'*P*dX;
    %     V3_dot(i)=[X(3:4,i);U-C_d*X(3:4,i)]'*P*dX+dX'*P*[X(3:4,i);U-C_d*X(3:4,i)];
    
    
    
    %Lyapunov function 4
    if flag
        if norm_dr>u_max/kr
            V4(i)=u_max*norm_dr+norm_dv^2/2-0.5*u_max^2/kr;
            V4_dot(i)=u_max*dr'*dv/norm_dr+dv_dot'*dv;
        else
            V4(i)=kr*norm_dr^2/2+norm_dv^2/2;
            V4_dot(i)=kr*dr'*dv+dv_dot'*dv;
        end
    else
        if norm_dr>(u_max1/kr)^(1/alpha_r)
            V4(i)=u_max1*norm_dr+norm_dv^2/2-alpha_r/(1+alpha_r)*u_max1^((1+alpha_r)/(alpha_r))/kr^(1/alpha_r);
            V4_dot(i)=u_max1*dr'*dv/norm_dr+dv_dot'*dv;
        else
            V4(i)=kr/(1+alpha_r)*norm_dr^(1+alpha_r)+norm_dv^2/2;
            V4_dot(i)=kr*dr'*dv*norm_dr^(alpha_r-1)+dv_dot'*dv;
        end
    end
    temp1=-kv*norm_dv^(alpha_v+1)*u_max2/norm_h2;
    %temp2=C_d*(u_max2-norm_h2)*dv'*X(3:4,i)/norm_h2;
    temp2=C_d*norm_v*(u_max2-norm_h2)*dv'*X(3:4,i)/norm_h2;
    %V4_dot0=dv'*((min(norm_h2,u_max2)*h2/norm_h2)-C_d*X(3:4,i))
    %V4_dot(i)=-C_d*norm_v^2;
    % V4_dot(i)=dvr*(norm_v^2+U'*dr);
    
    %Lyapunov function 5
    V5(i)=norm_dr;
    V5_dot(i)=dvr/norm_dr;
    V5_ddot(i)=(norm_dv^2+dv_dot'*dr)/norm_dr-dvr^2/norm_dr^3;
    
    t(i+1)=t(i)+dt;
    vel_norm(i+1)=norm(X(3:4,i+1));
    norm_X(i)=norm(X(:,i));
    norm_dX(i)=norm(X(:,i)-[xi;eta]);
    
    
    if norm(X(1:2,i)-xi)<1e-5
        break;
    end
end
%V3_dot = diff(V3);
% V3_dot(i+1) = V3_dot(i);

V1(i+1)=V1(i);
V1_dot(i+1)=V1_dot(i);
V1_ddot(i+1)=V1_ddot(i);

V2(i+1)=V2(i);
V2_dot(i+1)=V2_dot(i);
V2_ddot(i+1)=V2_ddot(i);

V3(i+1)=V3(i);
V3_dot(i+1)=V3_dot(i);

V4(i+1)=V4(i);
V4_dot(i+1)=V4_dot(i);

V5(i+1)=V5(i);
V5_dot(i+1)=V5_dot(i);
V5_ddot(i+1)=V5_ddot(i);

dist_norm(i+1)=dist_norm(i);
flag_sat(i+1)=flag_sat(i);
norm_U(i+1)=norm_U(i);
norm_X(i+1)=norm_X(i);
norm_dX(i+1)=norm_dX(i);

figure
plot(t,X(1,:),t,X(2,:))
legend('x','y')

figure
plot(t,vel_norm)
hold on;
plot(t,ones(size(t))*v_max,'r--')
plot(t,flag_sat,'g--')
ylabel('$||v||$')

figure
plot(t,V1,'b')
hold on
plot(t,V1_dot,'r')
plot(t,V1_ddot,'k')
ylabel('$V_1$')

figure
plot(t,V2,'b')
hold on
plot(t,V2_dot,'r')
plot(t,max(V2)*flag_sat,'g--')
plot(t,V2_ddot,'k')
ylabel('$V_2$')

figure
plot(t,V3,'b')
hold on
plot(t,V3_dot,'r')
plot(t,cm*flag_sat,'g--')
plot(t,cm*ones(size(t)),'m--')
%plot(t,c0*ones(size(t)),'b.-')
ylabel('$V_3$')

figure
plot(t,V4,'b')
hold on
plot(t,V4_dot,'r')
plot(t,max(V4)*flag_sat,'g--')
ylabel('$V_4$')

figure
plot(t,V5,'b')
hold on
plot(t,V5_dot,'r')
plot(t,V5_ddot,'k')
plot(t,R_max*ones(size(t)),'--')
plot(t,flag_sat,'g--')
ylabel('$V_5$')

figure
plot(t,norm_U)
hold on
plot(t,flag_sat,'g--')
ylabel('$||U||$')

figure
plot(t,norm_dX)
hold on
plot(t,flag_sat,'g--')
ylabel('$||dX||$')

if(1)
    figure
    hold on;
    
    %plot the circle outside which control is saturated
    plot(xi(1)+R_max*cos(0:pi/100:2*pi), xi(2)+R_max*sin(0:pi/100:2*pi),'r--')
    %
    plot(xi(1)+R_max_inv1*cos(0:pi/100:2*pi), xi(2)+R_max_inv1*sin(0:pi/100:2*pi),'m--')
    plot(xi(1)+R_max_inv2*cos(0:pi/100:2*pi), xi(2)+R_max_inv2*sin(0:pi/100:2*pi),'k--')
    xlabel('x');
    ylabel('y');
    %For Invariant level set of V=x'Px
    
    ind3=find(flag_sat==3);
    ind2=find(flag_sat==2);
    ind1=find(flag_sat==1);
    ind0=find(flag_sat==0);
    plot(X(1,ind3),X(2,ind3),'co','markersize',2);
    plot(X(1,ind2),X(2,ind2),'mo','markersize',2);
    plot(X(1,ind1),X(2,ind1),'ro','markersize',2);
    plot(X(1,ind0),X(2,ind0),'bo','markersize',2);
    if flag_sat(1)==1
        hand=plot(X(1,1),X(2,1),'rsquare','markersize',6);
        plot(X(1,1),X(2,1),'rsquare','markersize',8)
    else
        hand=plot(X(1,1),X(2,1),'bsquare','markersize',6);
        plot(X(1,1),X(2,1),'bsquare','markersize',8)
    end
    for j=2:10:i
        if flag_sat(i)==1
            set(hand,'XData',X(1,j),'YData',X(2,j),'color',[1,0,0])
        else
            set(hand,'XData',X(1,j),'YData',X(2,j),'color',[0,0,1])
        end
        drawnow;
    end
end
%%
syms x y
f=@(x,y) p1*x^2+p2*y^2+2*p12*x*y-c;

x0=-R_max:2*R_max/200:R_max;
y0=-v_max:2*v_max/200:v_max;
[x1,y1]=meshgrid(x0,y0);
v=p1*x1.^2+p2*y1.^2+2*p12*x1.*y1;
figure
fcontour(f)
max_V = -100;
for i =1:length(x0)
    for j=1:length(y0)
        x=x0(i);
        y=y0(i);
        V = p1*x^2+p2*y^2+2*p12*x*y;
        if V>max_V
            max_V = V;
        end
    end
end
%%
%Bound on time derivative of lyapunov function V=x'Px
v=0:v_max/100:v_max;
for i=1:length(v)
    V_dot_bound(i)=-v(i)^2 +((2*p1-2*p12*C_d)*sqrt(c/p1)-2*p2*u_max)*(c-(p1*R_max^2+p2*v(i)^2))/(2*p12*R_max)-2*p12*u_max*sqrt(c/p1);
end

figure
plot(v,V_dot_bound)


r0=R_max:(R_maxV-R_max)/100:R_maxV;
v0=0:v_max/100:v_max;
[r,v]=meshgrid(r0,v0);
V_dot=zeros(size(r));

for i=1:length(r)
    for j=1:length(v)
        if p1*r(i,j)^2+p2*v(i,j)^2+2*p12*r(i,j)*v(i,j)<c/1.5
            V_dot(i,j)=-v(i,j)^2 +((2*p1-2*p12*C_d)*r(i,j)-2*p2*u_max)*v(i,j)-2*p12*u_max*r(i,j);
        end
    end
end

figure
surf(r,v,V_dot)

%%
%plot y=(2*p2*x+Q0*x^2)/(2*p1*x-2*p12*x-2*p12*u_max)
x10=0:10/100:10; %v
x20=x10;  %eta
[x1,x2]=meshgrid(x10,x20);
%z=zeros(size(x1));
max_z=-inf;
for i=1:length(x10)
    for j=1:length(x20)
        if x1(i,j)<x2(i,j)
        z1(i,j)=-kv*(x2(i,j)-x1(i,j))^(alpha_v+1)*u_max2;
        z2(i,j)=C_d*x1(i,j)^2.*(u_max2-C_d*x1(i,j)^2-kv*(x1(i,j)+x2(i,j))^alpha_v)*(x1(i,j)-x2(i,j));
        z(i,j)=z1(i,j)+z2(i,j);
        if (z(i,j))>max_z
            max_z=z(i,j);
        end
        if max_z>0
            1
        end
        end
    end
end
figure
surf(x1,x2,z)
xlabel('v')
ylabel('$\eta$')

%%
x10=0:10/100:v_max/3; %eta
x20=x10;  %v_tilde
%v_max=20;
[x1,x2]=meshgrid(x10,x20);
%z=zeros(size(x1));
max_z=-inf;
for i=1:length(x10)
    for j=1:length(x20)
        if x2(i,j)<max(x1(i,j),v_max-x1(i,j)) && u_max2<abs(C_d*(x1(i,j)-x2(i,j))^2-kv*x2(i,j)^alpha_v)
        z1(i,j)=-kv*(x2(i,j))^(alpha_v+1)*u_max2;
        z2(i,j)=C_d*(x1(i,j)+x2(i,j))^4*x2(i,j)+C_d*kv*(x1(i,j)+x2(i,j))^2*x2(i,j)^(alpha_v+1)-u_max2*C_d*(x1(i,j)-x2(i,j))^2*x2(i,j);
       % z2(i,j)=z2(i,j)*x1(i,j)/(x1(i,j)+x2(i,j));
        z(i,j)=z1(i,j)+z2(i,j);
        if (z(i,j))>max_z
            max_z=z(i,j);
        end
        if max_z>0
            1
        end
        end
    end
end
figure
surf(x1,x2,z)
xlabel('v')
ylabel('$\eta$')