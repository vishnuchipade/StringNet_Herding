% ------------------------------
% Author: Vishnu S. Chipade
% Email: vishnuc@umich.edu
% ------------------------------

%SetPlotDefaults;

fontSize=18;
colors={[0,0,1],[1,0,0],[0,0,0],[0.5,0.5,0.5],[0,0.5,1],[1,0.5,1],[0.5,0.5,1],[1,.1,0.5],[0.1,0.9,0.4],[0.9,0.4,0.5],[0,0.5,0.4]};

circCos=cos(0:pi/100:2*pi);
circSin=sin(0:pi/100:2*pi);
%draw two lines intersecting each other
L1=[5,5;8,10]';
L2=[7.0,5;5,8]';
rho_D=.8;
theta1=atan2(L1(2,2)-L1(2,1),L1(1,2)-L1(1,1));
mL1=tan(theta1);
cL1=L1(2,1)-mL1*L1(1,1);
theta2=atan2(L2(2,2)-L2(2,1),L2(1,2)-L2(1,1));
mL2=tan(theta2);
cL2=L2(2,1)-mL2*L2(1,1);
x_int=(cL2-cL1)/(mL1-mL2);
y_int=mL1*x_int+cL1;
dtheta=abs(theta1-theta2);
intS=rho_D*sqrt(2/(1-abs(cos(dtheta))))
%Line 1 intersection points
dx=intS*abs(cos(theta1));
x11=x_int-dx;
y11=mL1*x11+cL1;
x12=x_int+dx;
y12=mL1*x12+cL1;
rec1=[x_int*ones(1,5);y_int*ones(1,5)]+[cos(theta1),-sin(theta1);sin(theta1),cos(theta1)]*([-intS,intS,intS,-intS,-intS;-rho_D,-rho_D,rho_D,rho_D,-rho_D])


%Line 2 intersection points
dx=intS*abs(cos(theta2));
x21=x_int-dx;
y21=mL2*x21+cL2;
x22=x_int+dx;
y22=mL2*x22+cL2;
rec2=[x_int*ones(1,5);y_int*ones(1,5)]+[cos(theta2),-sin(theta2);sin(theta2),cos(theta2)]*([-intS,intS,intS,-intS,-intS;-rho_D,-rho_D,rho_D,rho_D,-rho_D])


figure
hold all
plot(L1(1,:),L1(2,:),'color',colors{1})
text(L1(1,1)+0.1,L1(2,1),'$P_1$','fontsize',fontSize)
plot(L2(1,:),L2(2,:),'color',colors{1})
text(L2(1,1)-0.35,L2(2,1),'$P_2$','fontsize',fontSize)
%plot the intersecting portion
plot([x11,x12],[y11,y12],'color',colors{2})  %centerline
plot(x11+rho_D*circCos,y11+rho_D*circSin,'--','color',colors{6})  %end circle
plot(x12+rho_D*circCos,y12+rho_D*circSin,'--','color',colors{6})
plot(rec1(1,:),rec1(2,:),'--','color',colors{6})  %swapping rectangle


plot([x21,x22],[y21,y22],'color',colors{2})  %centerline
plot(x21+rho_D*circCos,y21+rho_D*circSin,'--','color',colors{6}) %end circle
plot(x22+rho_D*circCos,y22+rho_D*circSin,'--','color',colors{6})
plot(rec2(1,:),rec2(2,:),'--','color',colors{6})  %swapping rectangle

%Angle arc between two lines
plot(x_int+0.4*rho_D*cos(theta1:dtheta/20:theta2),y_int+0.4*rho_D*sin(theta1:dtheta/20:theta2),'color',colors{1})
%Text
text(x_int-0.1,y_int+.5,'$\theta$','fontsize',fontSize)  %theta
quiver(x_int+dx/7,mL1*(x_int+dx/7)+cL1,-cos(theta2),-sin(theta2),1.6,'color',colors{2}) 
text(6.82,5.82,'$\bar{P}$','color',colors{2},'fontsize',fontSize)
text(5.5,9.5,'$\bar{P}=\rho_d\sqrt{\frac{2}{1-\cos(\theta)}}$','fontsize',fontSize)
