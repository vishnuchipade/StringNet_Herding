%this function plots the beta agent position and velocity on superelliptic
% ------------------------------
% Author: Vishnu S. Chipade
% Email: vishnuc@umich.edu
% ------------------------------

%curves around rectangular obstacles

options = optimset('Display','off','MaxIter',1000);

rI=[0,-10]';
%rS=[-3,39]';
%rO=[6,21;8,-7]';
rS=[5,7]';
rO=[0,0;8,-7]';
ROS=norm(rS-rO(:,1));
NO=size(rO,2);
rotSense=[2,1,2];
RO=2;
w=sqrt(2)*RO;
w=4;
h=w*0.75;
RO_max=6;
k=1;
rA=[4,1]';
vA=[-2,4]';
rA0=[5,-1]';
vA0=[-3,4]';
clear  X Y
xl=10;
[X,Y]=meshgrid(-xl+rO(1,1):.5:xl+rO(1,2),-xl+rO(2,2):.5:xl+rO(2,1));

alpha=1;

Ktol=0.05;
Rmax=50;
rho_F=0.3; %formation radius
f1=@(x) x-0.5*(2*(w/2+rho_F)/w)^(2/(1-exp(-x)))-0.5*(2*(h/2+rho_F)/h)^(2/(1-exp(-x)))+1;
K0=fsolve(f1,rho_F+sqrt(w^2+h^2)/2,options);
%K0=0.8;
nmin=1/(1-exp(-alpha*K0));
n=nmin;
a=w/2*(2^(1/(2*n)));
b=h/2*(2^(1/(2*n)));
Ktol=K0;
EmAO=K0;
%Ktol=1.83;
pt1=27;
pt2=15;
Kmax=3;
Kmax2=abs(rA0(1)-rO(1,k))^(2*n)/a^(2*n)+abs(rA0(2)-rO(2,k))^(2*n)/b^(2*n)-1;
EuAO=Kmax2;
dKcube=(Kmax2-Kmax)^3;
A=2/dKcube;
B=-3*(Kmax2+Kmax)/dKcube;
C=6*Kmax2*Kmax/dKcube;
D=Kmax2^2*(Kmax2-3*Kmax)/dKcube;
%Point along the radial direction as the initial condition (rA should be in first quadrant of rO)
beta=atan2(rA(2)-rO(2),rA(1)-rO(1));
rB=((EmAO+1)/(abs(cos(beta))^(2*n)/a^(2*n)+abs(sin(beta))^(2*n)/b^(2*n)))^(1/(2*n));
rP0=[rB*cos(beta);rB*sin(beta)];
%Find the projection of rA on EmAO
f=@(r) [(r(1)-rO(1))^(2*n)/a^(2*n)+(r(2)-rO(2))^(2*n)/b^(2*n)-EmAO-1; a^(2*n)*(r(2)-rO(2))^(2*n-1)*(rA(1)-r(1))-b^(2*n)*(r(1)-rO(1))^(2*n-1)*(rA(2)-r(2))]
f=@(r) [abs(r(1)-rO(1,k))^(2*n)/a^(2*n)+abs(r(2)-rO(2,k))^(2*n)/b^(2*n)-EmAO-1; a^(2*n)*sign(r(2)-rO(2,k))*abs(r(2)-rO(2,k))^(2*n-1)*(rA(1)-r(1))-b^(2*n)*sign(r(1)-rO(1,k))*abs(r(1)-rO(1,k))^(2*n-1)*(rA(2)-r(2))];
rAProj=fsolve(f,rP0,options)
%Tangent at the projection point
beta_barP=atan2(b^(2*n)*(rAProj(1)-rO(1))^(2*n-1),-a^(2*n)*(rAProj(2)-rO(2))^(2*n-1));
rTP=[cos(beta_barP);sin(beta_barP)];
vAProj=rTP'*vA*rTP;


%Find the straight line portion of the obstacle boundary
k=1;
BetaAv0=atan2(vA0(2),vA0(1)); %Orientation while entering the obstacle field
SpeedA0=norm(vA0); %Speed while entering the obstacle field
mAv0=tan(BetaAv0);
beta=atan(sign(-b^(2*n)/(mAv0*a^(2*n)))*abs(-b^(2*n)/(mAv0*a^(2*n)))^(1/(2*n-1)));
rB=((EmAO+1)/(abs(cos(beta))^(2*n)/a^(2*n)+abs(sin(beta))^(2*n)/b^(2*n)))^(1/(2*n));
rP0=rO(1:2,k)+[rB*cos(beta);rB*sin(beta)];
%rAProj=rP0;
%beta_bar=atan2(b^(2*n)*sign(rAProj(1,k)-rO(1,k))*abs(rAProj(1,k)-rO(1,k))^(2*n-1),-a^(2*n)*sign(rAProj(2,k)-rO(2,k))*abs(rAProj(2,k)-rO(2,k))^(2*n-1));
crossProd=cross([cos(BetaAv0);sin(BetaAv0);0],[rO(1:2,k)-rA0;0]);
if crossProd(3)>0
    crossProd2=cross([rP0-rO(1:2,k);0],[rO(1:2,k)-rA0;0]);
    if crossProd2(3)<0
        rP0=rO(1:2,k)-[rB*cos(beta);rB*sin(beta)];
    end
else
    crossProd2=cross([rP0-rO(1:2,k);0],[rO(1:2,k)-rA0;0]);
    if crossProd2(3)>0
        rP0=rO(1:2,k)-[rB*cos(beta);rB*sin(beta)];
    end
end
cAv0=rP0(2)-mAv0*rP0(1);


colors={[0.9 0.1 0.1],[0.1 0.7 0.1],[0.1 0.7 0.1]}

figure
hold all;
xsup=a*(EmAO+1)^(1/(2*n))*cos(0:pi/50:pi/2).^(1/n);
ysup=b*(EmAO+1)^(1/(2*n))*sin(0:pi/50:pi/2).^(1/n);
xsup=[xsup,-flip(xsup),-xsup,flip(xsup)];
ysup=[ysup,flip(ysup),-ysup,-flip(ysup)];
xsup=xsup+rO(1,k);
ysup=ysup+rO(2,k);
plot(xsup,ysup,'color',[0.5,0,1],'linewidth',0.5)
%set(ez,'color',[0.5,0,1])
if(0)
    xsup=a*(EbarAO+1)^(1/(2*n))*cos(0:pi/50:pi/2).^(1/n);
    ysup=b*(EbarAO+1)^(1/(2*n))*sin(0:pi/50:pi/2).^(1/n);
    xsup=[xsup,-flip(xsup),-xsup,flip(xsup)];
    ysup=[ysup,flip(ysup),-ysup,-flip(ysup)];
    xsup=xsup+rO(1,k);
    ysup=ysup+rO(2,k);
    plot(xsup,ysup,'color',[0.3,.91,0],'linewidth',0.5)
end
xsup=a*(EuAO+1)^(1/(2*n))*cos(0:pi/50:pi/2).^(1/n);
ysup=b*(EuAO+1)^(1/(2*n))*sin(0:pi/50:pi/2).^(1/n);
xsup=[xsup,-flip(xsup),-xsup,flip(xsup)];
ysup=[ysup,flip(ysup),-ysup,-flip(ysup)];
xsup=xsup+rO(1,k);
ysup=ysup+rO(2,k);
plot(xsup,ysup,'color',[0,.0,0],'linewidth',0.5)
%position of A
plot(rA(1),rA(2),'o','color',colors{1},'markersize',5);
%Projection
plot(rAProj(1),rAProj(2),'square','color',colors{2});
legend(['$E_o=\xi_o^m=$',num2str(EmAO,3)],['$E_o=\xi_o^u=$',num2str(EuAO,3)],['$\mathbf{r}_{a}$'],['$\mathbf{r}_{\beta}$'],'location','northwest','AutoUpdate','off' )

rectangle('position',[rO(1,1)-w(k)/2,rO(2,1)-h(k)/2,w(k),h(k)],'FaceColor',[0.5 .5 .5])
hold on;
viscircles([rO(1,1)-w(k)/2,rO(2,1)-h(k)/2;rO(1,1)+w(k)/2,rO(2,1)-h(k)/2;rO(1,1)+w(k)/2,rO(2,1)+h(k)/2;rO(1,1)-w(k)/2,rO(2,1)+h(k)/2],rho_F*ones(1,4),'color','r','LineStyle','--','linewidth',.5)
hold on;
line([rO(1,1)-w(k)/2-rho_F,rO(1,1)+w(k)/2+rho_F],[rO(2,1)-h(k)/2-rho_F,rO(2,1)-h(k)/2-rho_F],'color','r','linestyle','--','linewidth',.5)
hold on
line([rO(1,1)-w(k)/2-rho_F,rO(1,1)-w(k)/2-rho_F],[rO(2,1)-h(k)/2-rho_F,rO(2,1)+h(k)/2+rho_F],'color','r','linestyle','--','linewidth',.5)
hold on;
line([rO(1,1)-w(k)/2-rho_F,rO(1,1)+w(k)/2+rho_F],[rO(2,1)+h(k)/2+rho_F,rO(2,1)+h(k)/2+rho_F],'color','r','linestyle','--','linewidth',.5)
hold on
line([rO(1,1)+w(k)/2+rho_F,rO(1,1)+w(k)/2+rho_F],[rO(2,1)-h(k)/2-rho_F,rO(2,1)+h(k)/2+rho_F],'color','r','linestyle','--','linewidth',.5)
hold on;

headSize=2;
fontSize=15;
%reference line
plot([rO(1),rO(1)+5],[rO(2),rO(2)],'b--')
%position of A
plot(rA(1),rA(2),'o','color',colors{1},'markersize',5);
if(0)
plot(rA0(1),rA0(2),'o','color',colors{1},'markersize',5);
end
%line to rA
plot([rO(1),rA(1)],[rO(2),rA(2)],'color',colors{1})
%velocity of A
quiver(rA(1),rA(2),vA(1),vA(2),0.25,'linewidth',2,'color',colors{1},'MaxHeadSize',headSize)
text(rA(1)-0.2,rA(2)+.5,'$\mathbf{v}_a$','fontsize',fontSize)
if(0)
%initial vA
quiver(rA0(1),rA0(2),vA0(1),vA0(2),0.25,'linewidth',2,'color',colors{1},'MaxHeadSize',headSize)
end
%Projection
plot(rAProj(1),rAProj(2),'square','color',colors{2});
%line to rAProj
plot([rO(1),rAProj(1)],[rO(2),rAProj(2)],'color',colors{2})
%Projected velocity
quiver(rAProj(1),rAProj(2),vAProj(1),vAProj(2),0.25,'linewidth',2,'color',colors{2},'MaxHeadSize',headSize)
text(rAProj(1),rAProj(2)+.5,'$\mathbf{v}_{\beta}$','fontsize',fontSize)

%Plot the line tangent to superelliptic curve and parallel to initial
%velocity
if(0)
xL=[1:0.01:rP0(1)];
yL=mAv0.*xL+cAv0;
hold on;
plot(xL,yL,'color',[0.5,0,1]);
end

