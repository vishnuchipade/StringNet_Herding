set(groot,'defaulttextinterpreter','latex');
% ------------------------------
% Author: Vishnu S. Chipade
% Email: vishnuc@umich.edu
% ------------------------------

set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultTextInterpreter','latex');
set(groot, 'DefaultLineLineWidth', 2);
set(groot, 'DefaultAxesFontSize', 18);

colors={[0,0,1],[1,0,0],[0,0,0],[0.5,0.5,0.5],[0,0.5,1],[1,0.5,1],[0.5,0.5,1],[1,.1,0.5],[0.1,0.9,0.4],[0.9,0.4,0.5],[0,0.5,0.4],[1,0.4,0.4]};
fontSize=20;
markSize=4;

circCos=cos(0:pi/50:2*pi);
circSin=sin(0:pi/50:2*pi);

rho_ac=8;
b_d=.5;
rho_ac_prime=rho_ac+b_d;
ND=6;
rho_sn=rho_ac_prime/cos(pi/ND);

figure('units','normalized','outerposition',[.1 0.1 .7 .8]);
hold all
plot(rho_ac*circCos,rho_ac*circSin,'r-')
plot(rho_ac_prime*circCos,rho_ac_prime*circSin,'r--')
plot(rho_sn*circCos,rho_sn*circSin,'b-')
axis equal;

%Draw stringNet
dr=[0.4,0;0.4,0.4;-1.5,0.5;-1.5,-0.3;-0.3,-1.15;0.2,-.75]';
for i=1:ND
    theta=(i-1)*2*pi/ND;
    theta2=(i)*2*pi/ND;
    rD(:,i)=[rho_sn*cos(theta),rho_sn*sin(theta)]';
    plot(rD(1,i),rD(2,i),'bo','markersize',markSize)
    text(rho_sn*cos(theta)+dr(1,i),rho_sn*sin(theta)+dr(2,i),['$\mathbf{\xi}_{',num2str(i),'}^e$'],'fontsize',fontSize)
    plot(rho_sn*[cos(theta),cos(theta2)],rho_sn*[sin(theta),sin(theta2)],'b-')
end

plot([0,rD(1,2)],[0,rD(2,2)],'b-')
text(rD(1,2)/2,rD(1,2)/2+1,'$\rho_{sn}$','fontsize',fontSize,'color','b')

theta=.9*pi;
plot([0,rho_ac*cos(theta)],[0,rho_ac*sin(theta)],'r-')
text(rho_ac/2*cos(theta)-.3,rho_ac/2*sin(theta)-1,'$\rho_{ac}$','fontsize',fontSize,'color','r')

theta=atan2(0.5*(rD(2,3)+rD(2,2)),0.5*(rD(1,3)+rD(1,2)));
plot([0,rho_ac_prime*cos(theta)],[0,rho_ac_prime*sin(theta)],'r--')
text(rho_ac_prime/2*cos(theta)-5,rho_ac_prime/2*sin(theta)-.2,'$\rho_{ac}+b_d$','fontsiz',fontSize,'color','r')
dp=1.1;
rPoint(:,1)=(rho_ac_prime-dp)*[cos(theta),sin(theta)]';
rPoint(:,2)=rPoint(:,1)+dp*[cos(theta-pi/2),sin(theta-pi/2)]';
rPoint(:,3)=rPoint(:,2)+dp*[cos(theta),sin(theta)]';
plot(rPoint(1,:),rPoint(2,:),'r--','linewidth',1.3)

text(-3,-5,'$\rho_{sn}=\rho_{sn}^{circum}$','fontsize',fontSize,'color','b')


text(-5,12.5,'$||\mathbf{\xi}_{2}^e-\mathbf{\xi}_{3}^e||=R_d^{d,s}-2b_d$','fontsize',fontSize,'color','b')
quiver(0,11.5,0,-2,1.35,'linewidth',1.5,'maxheadsize',.95,'color','b')

%Plot for the other case with rho_sn=rho_sn^max-b_d
rC=[25,0]';
plot(rC(1)+rho_ac*circCos,rC(2)+rho_ac*circSin,'r-')
plot(rC(1)+rho_ac_prime*circCos,rC(2)+rho_ac_prime*circSin,'r--')
plot(rC(1)++rho_sn*circCos,rC(2)+rho_sn*circSin,'b-')
axis equal;

%Draw stringNet
dr=[0.4,0;0.4,0.4;-1.5,0.5;-1.5,-0.3;-0.3,-1.15;0.2,-.75]';
for i=1:ND
    theta=(i-1)*2*pi/ND;
    theta2=(i)*2*pi/ND;
    rD(:,i)=[rho_sn*cos(theta),rho_sn*sin(theta)]';
    plot(rC(1)+rD(1,i),rC(2)+rD(2,i),'bo','markersize',markSize)
    text(rC(1)+rho_sn*cos(theta)+dr(1,i),rC(2)+rho_sn*sin(theta)+dr(2,i),['$\mathbf{\xi}_{',num2str(i),'}^e$'],'fontsize',fontSize)
    plot(rC(1)+rho_sn*[cos(theta),cos(theta2)],rC(2)+rho_sn*[sin(theta),sin(theta2)],'b-')
end

plot(rC(1)+[0,rD(1,2)],rC(2)+[0,rD(2,2)],'b-')
text(rC(1)+rD(1,2)/2,rC(2)+rD(1,2)/2+1,'$\rho_{sn}$','fontsize',fontSize,'color','b')

theta=.9*pi;
plot(rC(1)+[0,rho_ac*cos(theta)],rC(2)+[0,rho_ac*sin(theta)],'r-')
text(rC(1)+rho_ac/2*cos(theta)-.3,rC(2)+rho_ac/2*sin(theta)-1,'$\rho_{ac}$','fontsize',fontSize,'color','r')

theta=atan2(0.5*(rD(2,3)+rD(2,2)),0.5*(rD(1,3)+rD(1,2)));
plot(rC(1)+[0,rho_ac_prime*cos(theta)],rC(2)+[0,rho_ac_prime*sin(theta)],'r--')
text(rC(1)+rho_ac_prime/2*cos(theta)-5,rC(2)+rho_ac_prime/2*sin(theta)-.2,'$\rho_{ac}+b_d$','fontsiz',fontSize,'color','r')
dp=1.1;
rPoint(:,1)=(rho_ac_prime-dp)*[cos(theta),sin(theta)]';
rPoint(:,2)=rPoint(:,1)+dp*[cos(theta-pi/2),sin(theta-pi/2)]';
rPoint(:,3)=rPoint(:,2)+dp*[cos(theta),sin(theta)]';
plot(rC(1)+rPoint(1,:),rC(2)+rPoint(2,:),'r--','linewidth',1.3)

text(rC(1)+-4.5,rC(2)+-5,'$\rho_{sn}=\bar{\rho}_{sn}-b_d$','fontsize',fontSize,'color','b')

text(rC(1)-5,rC(2)+12.5,'$||\mathbf{\xi}_{2}^e-\mathbf{\xi}_{3}^e||<R_d^{d,s}-2b_d$','fontsize',fontSize,'color','b')
quiver(rC(1)+0,rC(2)+11.5,0,-2,1.35,'linewidth',1.5,'maxheadsize',.95,'color','b')

text(-5,-12,'$(i)\; \rho_{sn}^{circum}<\bar{\rho}_{sn}-b_d $','fontsize',fontSize)
text(20,-12,'$(ii) \;\rho_{sn}^{circum}\ge\bar{\rho}_{sn}-b_d $','fontsize',fontSize)

xlim([-13,36])

