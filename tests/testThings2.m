% ------------------------------
% Author: Vishnu S. Chipade
% Email: vishnuc@umich.edu
% ------------------------------
 
global R s2 vmc


R=15;
qm=1;
vm=4;
vmc=sqrt(qm*R);
s1=.1*R; s3=.1*R;  %less than R/2
s2=4*R;
v1m=sqrt(2*qm*s1);
v2m=sqrt(2*qm*s3);
v1=0:v1m/100:v1m;
v2=0:v2m/100:v2m;

v_half=@(x) vmc*sqrt(tanh(s2/R+0.5*atanh(x(1)^2/vmc^2)+0.5*atanh(x(2)^2/vmc^2)))
f=@(x) 1/(2*k2)*(log((v_half(x)+vmc)/(vmc-v_half(x)))+2*atan(v_half(x)/vmc))-1/(4*k2)*(log((x(1)+vmc)/(vmc-x(1)))+2*atan(x(1)/vmc))-1/(4*k2)*(log((x(2)+vmc)/(vmc-x(2)))+2*atan(x(2)/vmc))

Aeq=[1,0,;0,0];
beq=[2.3;0];
 nonlcon=@conc;
 v0=fmincon(f,[0.1,0.1],[],[],Aeq,beq,[0,0],[vmc,vmc],nonlcon)
 

 
%%
%Check intersection of circular robot on a circle and a line
rVC=[-100,-60]';
mL1=1;
cL1=5;
x1=-1.0;
x2=0.0;
r2=[x1,mL1*x1+cL1]';
r1=[x2,mL1*x2+cL1]';
r2=[-93.14,-57.80]';
r1=[0,-350]';
dr=r2-r1;
L1=norm(dr);
drx=dr(1);
dry=dr(2);
dr_hat=dr/norm(dr);
rho_safe=7.2;%2.5*sqrt(2);
rho_D=0.5;

segType=1;
mL1=dry/drx;
cL1=r1(2)-mL1*r1(1);
ang1=-7*pi/180;
ang2=45*pi/180;
r21=rVC+rho_safe*[cos(ang1),sin(ang1)]';
r22=rVC+rho_safe*[cos(ang2),sin(ang2)]';
% r21=[-92.8,-60]';
% r22=[-92.87,-59.07]';
interSec=interSecLineCirc(r1,r2,mL1,cL1,L1,dr_hat,drx,dry,r21,r22,rVC,segType);
plotIntersectionLineCirc(r1,r2,r21,r22,rVC,segType,interSec)

%%
if rc(2)-mL1*rc(1)-cL1>0   %if the center is above the path line
    cL11=cL1+2*rho_D*sqrt(1+mL1^2);
else
    cL11=cL1-2*rho_D*sqrt(1+mL1^2);
end
circCos=cos(0:pi/100:2*pi);
circSin=sin(0:pi/100:2*pi);
r1_p(1,1)=(r1(1)+mL1*(r1(2)-cL11))/(1+mL1^2);
r1_p(2,1)=r1(2)+(mL1*r1(1)+cL11-r1(2))/(1+mL1^2);
[lambda]=lambdaInterSecCircLine(rc,rhoc,r1_p,mL1,cL11,drx,dry);
x_int1=r1_p(1)+lambda(1)*drx; 
x_int2=r1_p(1)+lambda(2)*drx; 
y_int1=r1_p(2)+lambda(1)*dry;
y_int2=r1_p(2)+lambda(2)*dry;
%Projection on main path line segment
    x_intp1=(x_int1+mL1*(y_int1-cL1))/(1+mL1^2);
    y_intp1=y_int1+(mL1*x_int1+cL1-y_int1)/(1+mL1^2);
    x_intp2=(x_int2+mL1*(y_int2-cL1))/(1+mL1^2);
    y_intp2=y_int2+(mL1*x_int2+cL1-y_int2)/(1+mL1^2);
figure
plot(rc(1)+rhoc*circCos,rc(2)+rhoc*circSin)
hold on
plot([r1(1),r2(1)],[r1(2),r2(2)])
%plot intersection points and corresponding line
x0=[r1(1):(r2(1)-r1(1))/100:r2(1)];
y0=mL1*x0+cL11;
plot(x0,y0,r1_p(1),r1_p(2),'o')
plot(x_int1,y_int1,'bo',x_int2,y_int2,'bo')
%Plot line segment in conflict
plot([x_intp1,x_intp2],[y_intp1,y_intp2],'m')

%Plot circles
plot(x_int1+rho_D*circCos,y_int1+rho_D*circSin,'--');
plot(x_intp1+rho_D*circCos,y_intp1+rho_D*circSin,'--');
plot(x_int2+rho_D*circCos,y_int2+rho_D*circSin,'--');
plot(x_intp2+rho_D*circCos,y_intp2+rho_D*circSin,'--');



%%
%%Double integrator with drag term
X(:,1)=[1,.5,0,0]';
Cd=0.9;
K=12;
t(1)=0;
dt=0.001;
nSim=10000;

for i=1:nSim
    X(:,i+1)=X(:,i)+dt*[X(3:4,i); -K*X(1:2,i)-Cd*X(3:4,i)*norm(X(3:4,i))];
    t(i+1)=t(i)+dt;
end

figure
subplot(2,1,1)
plot(t,X(1,:),'b',t,X(2,:),'r')
ylabel('xs')
subplot(2,1,2)
plot(t,X(3,:),'b',t,X(4,:),'r')
xlabel('time')
ylabel('v')




%%
%Function definition
function [c,ceq]= conc(x)
global R s2 vmc
 c=R/2*abs(atanh(x(1)^2/vmc^2)-atanh(x(2)^2/vmc^2))-s2;
 ceq=[];
end


