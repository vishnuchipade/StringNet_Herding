% ------------------------------
% Author: Vishnu S. Chipade
% Email: vishnuc@umich.edu
% ------------------------------

%SetPlotDefaults;
set(groot,'defaulttextinterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

%%Draw defenders with string and attackers and the projection of the
%%attacker on the string

rD=[12,12;15,20;22,21;28,15;19,8]';
ND=size(rD,2);
vD=repmat([2,1]',1,ND);
rhoD=1;

rA=[19,15;20,18;24,16;22,13;18,11;]';
vA=[-3.5,2;-3.2,2;-3.,2;-2.5,2;-4.9,2;]';
NA=size(rA,2);
rhoA=rhoD;

circCos=cos(0:pi/100:2*pi);
circSin=sin(0:pi/100:2*pi);

colors={[0,0,1],[1,0,0.5],[0,0.5,0.5],[0.5,0.8,0.2]}

fontSize=18;
headSize=1;
markerSize=2;

%%
figure
hold all;
axis equal
%Initial string chain
if ND>1
    WDString=zeros(ND);
    NND=2;
    NAC=ND;
    temp=[1:NAC,1:NAC];
    for j=1:NAC
        WDString(j,temp([j+1:j+NND/2,j+NAC-NND+1:j+NAC-NND/2]))=ones(1,NND);
    end
else
    WDString=0;
end

%Plot defenders
for j=1:ND
    fill(rD(1,j)+rhoD*circCos,rD(2,j)+rhoD*circSin,'b');
    text(rD(1,j)-0.5,rD(2,j),['$\mathcal{D}_{',num2str(j),'}$'],'color',[1,1,1],'fontsize',fontSize);
     thetav=atan2(vD(2,j),vD(1,j));
    quiver(rD(1,j)+rhoA*cos(thetav),rD(2,j)+rhoA*sin(thetav),vD(1,j),vD(2,j),0.5,'color',[0,0,1],'MaxHeadSize',0.9*headSize);
end
%plot string lines joining the defenders
WDString_temp=WDString;
countAPS=0;
for j=1:ND
    for jj=find(WDString_temp(j,:)==1)
        theta=atan2(rD(2,jj)-rD(2,j),rD(1,jj)-rD(1,j));
        countAPS=countAPS+1;
        mDD=tan(theta);
        cDD=rD(2,j)-mDD*rD(1,j);
        rAProjS(1,countAPS)=(mDD*rA(2)+rA(1)-mDD*cDD)/(1+mDD^2);
        rAProjS(2,countAPS)=mDD*rAProjS(1,countAPS)+cDD; 
        thetaProjS(countAPS)=theta+pi/2;
                
        plot([rD(1,j)+rhoD*cos(theta),rD(1,jj)+rhoD*cos(theta+pi)],[rD(2,j)+rhoD*sin(theta),rD(2,jj)+rhoD*sin(theta+pi)],'b','linewidth',2)
        plot([rD(1,j)+rhoD*cos(theta),rD(1,jj)+rhoD*cos(theta+pi)],[rD(2,j)+rhoD*sin(theta),rD(2,jj)+rhoD*sin(theta+pi)],'w--','linewidth',1)
        text((rD(1,j)+rD(1,jj))/2,(rD(2,j)+rD(2,jj))/2,['$\mathcal{B}_{',num2str(j),'}$'],'fontsize',fontSize)
        %Remove the strings from plot
        WDString_temp(j,jj)=0;
        WDString_temp(jj,j)=0;
    end
end

%Plot the attackers
for i=1:NA
    fill(rA(1,i)+rhoA*circCos,rA(2,i)+rhoA*circSin,'r');
    text(rA(1,i)-0.5,rA(2,i),['$\mathcal{A}_{',num2str(i),'}$'],'color',[1,1,1],'fontsize',fontSize);
    thetav=atan2(vA(2,i),vA(1,i));
    quiver(rA(1,i)+rhoA*cos(thetav),rA(2,i)+rhoA*sin(thetav),vA(1,i),vA(2,i),.3,'color',[1,0,0],'MaxHeadSize',0.9*headSize);
end

%Normal projection of A on the strings joining the defenders 
 for i=1:1
     theta=thetaProjS(i);
    plot([rA(1,i)+rhoA*cos(theta),rAProjS(1,i)],[rA(2,i)+rhoA*sin(theta),rAProjS(2,i)],'r--') 
    plot(rAProjS(1,i),rAProjS(2,i),'ro')
    text(0.3*rA(1,i)+0.7*rAProjS(1,i),0.3*rA(2,i)+0.7*rAProjS(2,i)+.5,'$R_{a1}^{b1}$','color',[0,0,0],'fontsize',fontSize)
 end

 
 
 
%%
%%Draw defenders with string in open formation and closed formation
rhoD=.5;

rA=[16,16;20,18;22,21;22.1,14;21,11.5;]';
NA=size(rA,2);
vA=repmat([-3.5,-2.5]',1,NA);

rhoA=rhoD;
rAcm=sum(rA,2)/NA;
rho_Acon=7;
rho_sn=10;
rho_sn_max=13;

thetaA=atan2(vA(2),vA(1))+38*pi/180;
RD0=rho_sn;
b=.5;
dtheta=acos(1-(2*b)^2/2/rho_sn^2);
ND=5;
for j=1:ND
    thetaD=thetaA-pi/2+pi*(j-1)/(ND-1);
rD(:,j)=[6,5.6]'+RD0*[cos(thetaD);sin(thetaD)];
end
ND=size(rD,2);
vD=0*repmat([2,1]',1,ND);


for j=1:ND
    thetaD=18*pi/180+thetaA-28*pi/180-pi+2*pi*(j-1)/ND;%+(j)*dtheta/ND;
rD_des(:,j)=rAcm+RD0*[cos(thetaD);sin(thetaD)];
end

circCos=cos(0:pi/100:2*pi);
circSin=sin(0:pi/100:2*pi);



figure('units','normalized','outerposition',[.1 0.1 .42 .7])
hold all;
grid on
grid minor
axis equal
%set(gca,'visible','off')
%set(gca,'xtick',[],'ytick',[])
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];

%Initial string chain
if ND>1
    WDString=zeros(ND);
    NND=2;
    NAC=ND;
    temp=[1:NAC,1:NAC];
    for j=1:NAC
        WDString(j,temp([j+1:j+NND/2,j+NAC-NND+1:j+NAC-NND/2]))=ones(1,NND);
    end
else
    WDString=0;
end
WDString(1,ND)=0;
WDString(ND,1)=0;

%Plot the attackers
for i=1:NA
    plot(rA(1,i)+rhoA*circCos,rA(2,i)+rhoA*circSin,'color',colors{2});
    fill(rA(1,i)+rhoA*circCos,rA(2,i)+rhoA*circSin,colors{2});
    thetaA=atan2(vA(2,i),vA(1,i));
    rHead=repmat(rA(:,i),1,4)+0.8*rhoA*[cos(thetaA),sin(thetaA); -sin(thetaA),cos(thetaA); sin(thetaA),-cos(thetaA);cos(thetaA),sin(thetaA)]';
    plot(rHead(1,:),rHead(2,:),'color',colors{2})  
    fill(rHead(1,:),rHead(2,:),[1,1,1])  
    text(rA(1,i)+0.8,rA(2,i),['$\mathcal{A}_',num2str(i),'$'],'fontsize',fontSize,'color',colors{2})
end

%plot the CoM of attackers
plot(rAcm(1),rAcm(2),'rsquare')
text(rAcm(1)-1.9,rAcm(2)+1,'$\mathbf{r}_{ac}$','fontsize',fontSize)
plot(rAcm(1)+rho_Acon*circCos,rAcm(2)+rho_Acon*circSin,'r')
R=rho_Acon+b;
plot(rAcm(1)+R*circCos,rAcm(2)+R*circSin,'r--')

%Plot radius of the connectivity region
theta=thetaA+40*pi/180;
plot([rAcm(1),rAcm(1)+rho_Acon*cos(theta)],[rAcm(2),rAcm(2)+rho_Acon*sin(theta)],'r')
text([rAcm(1)+rAcm(1)+rho_Acon*cos(theta)]/2-2,[rAcm(2)+rAcm(2)+rho_Acon*sin(theta)]/2+0.2,'$\rho_{ac}$', 'fontsize',fontSize,'color',[1,0,0])
theta=0;
% plot([rAcm(1),rAcm(1)+R*cos(theta)],[rAcm(2),rAcm(2)+R*sin(theta)],'r--')
% text([rAcm(1)+rAcm(1)+R*cos(theta)]/2-2.5,[rAcm(2)+rAcm(2)+R*sin(theta)]/2+.8,'$\rho_{ac}+b_d$', 'fontsize',fontSize)

%Plot defenders in open formation
for j=1:ND
    fill(rD(1,j)+rhoD*circCos,rD(2,j)+rhoD*circSin,'b');
    %text(rD(1,j)-0.5,rD(2,j)-0.5,['$\mathcal{D}_{',num2str(j),'}$'],'color',[1,1,1],'fontsize',fontSize);
    %text(rD(1,j)-2.5/j,rD(2,j)-j*0.3,['${\xi}_{',num2str(j),'}^o$'],'color',[0,0,1],'fontsize',fontSize);
     thetav=atan2(vD(2,j),vD(1,j));
    quiver(rD(1,j)+rhoA*cos(thetav),rD(2,j)+rhoA*sin(thetav),vD(1,j),vD(2,j),'color',[0,0,1]);
end
 text(rD(1,1)-1.3,rD(2,1)+1.3,['${\xi}_{',num2str(1),'}^{g}$'],'color',[0,0,1],'fontsize',fontSize);
 text(rD(1,2)-1.2,rD(2,2)-1.6,['${\xi}_{',num2str(2),'}^{g}$'],'color',[0,0,1],'fontsize',fontSize);
 text(rD(1,3)-1.9,rD(2,3)-.6,['${\xi}_{',num2str(3),'}^{g}$'],'color',[0,0,1],'fontsize',fontSize);
 text(rD(1,4)-.5,rD(2,4)-1.3,['${\xi}_{',num2str(4),'}^{g}$'],'color',[0,0,1],'fontsize',fontSize);
 text(rD(1,5)+1,rD(2,5)-.3,['${\xi}_{',num2str(5),'}^{g}$'],'color',[0,0,1],'fontsize',fontSize);
%Plot defenders' desired positions
for j=1:ND
    fill(rD_des(1,j)+rhoD*circCos,rD_des(2,j)+rhoD*circSin,[0,0.5,0.5]);
   % text(rD_des(1,j)-0.5,rD_des(2,j),['$\mathcal{D}_{',num2str(j),'}$'],'color',[1,1,1],'fontsize',fontSize);
    %text(rD_des(1,j)-2.5/j,rD_des(2,j)-j*0.3,['${\xi}_{',num2str(j),'}^e$'],'color',[0,0.5,0.5],'fontsize',fontSize);
     thetav=atan2(vD(2,j),vD(1,j));
    quiver(rD(1,j)+rhoA*cos(thetav),rD(2,j)+rhoA*sin(thetav),vD(1,j),vD(2,j),'color',[0,0,1]);
end
text(rD_des(1,1)-.3,rD_des(2,1)+1.3,['${\xi}_{',num2str(1),'}^e$'],'color',[0,0.5,0.5],'fontsize',fontSize);
 text(rD_des(1,2)-1.7,rD_des(2,2)-.6,['${\xi}_{',num2str(2),'}^e$'],'color',[0,0.5,0.5],'fontsize',fontSize);
 text(rD_des(1,3)-1.9,rD_des(2,3)-.6,['${\xi}_{',num2str(3),'}^e$'],'color',[0,0.5,0.5],'fontsize',fontSize);
 text(rD_des(1,4)-.5,rD_des(2,4)-1.3,['${\xi}_{',num2str(4),'}^e$'],'color',[0,0.5,0.5],'fontsize',fontSize);
 text(rD_des(1,5)+.1,rD_des(2,5)-1.3,['${\xi}_{',num2str(5),'}^e$'],'color',[0,0.5,0.5],'fontsize',fontSize);
%plot string lines joining the defenders
WDString_temp=WDString;
countAPS=0;
for j=1:ND
    for jj=find(WDString_temp(j,:)==1)
        theta=atan2(rD(2,jj)-rD(2,j),rD(1,jj)-rD(1,j));
        countAPS=countAPS+1;
        mDD=tan(theta);
        cDD=rD(2,j)-mDD*rD(1,j);
        rAProjS(1,countAPS)=(mDD*rA(2)+rA(1)-mDD*cDD)/(1+mDD^2);
        rAProjS(2,countAPS)=mDD*rAProjS(1,countAPS)+cDD; 
        thetaProjS(countAPS)=theta+pi/2;
                
        plot([rD(1,j)+rhoD*cos(theta),rD(1,jj)+rhoD*cos(theta+pi)],[rD(2,j)+rhoD*sin(theta),rD(2,jj)+rhoD*sin(theta+pi)],'b','linewidth',2)
        plot([rD(1,j)+rhoD*cos(theta),rD(1,jj)+rhoD*cos(theta+pi)],[rD(2,j)+rhoD*sin(theta),rD(2,jj)+rhoD*sin(theta+pi)],'w--','linewidth',1)
        text((rD(1,j)+rD(1,jj))/2,(rD(2,j)+rD(2,jj))/2,['$\mathcal{B}_{',num2str(j),'}$'],'fontsize',fontSize)
        %Remove the strings from plot
        WDString_temp(j,jj)=0;
        WDString_temp(jj,j)=0;
    end
end
%Initial string chain
if ND>1
    WDString=zeros(ND);
    NND=2;
    NAC=ND;
    temp=[1:NAC,1:NAC];
    for j=1:NAC
        WDString(j,temp([j+1:j+NND/2,j+NAC-NND+1:j+NAC-NND/2]))=ones(1,NND);
    end
else
    WDString=0;
end


%Plot the center and orientation of the open formation
rDc=(rD(:,1)+rD(:,end))/2;
dc=0.6;
rDc1=dc*rD(:,1)+(1-dc)*rD(:,end);
rDc2=(1-dc)*rD(:,1)+dc*rD(:,end);
plot(rDc(1),rDc(2),'bo')
text(rDc(1),rDc(2)-1,'$\mathbf{r}_{df}$','fontsize',fontSize)
plot([rDc(1),rD(1,2)],[rDc(2),rD(2,2)],'b')
text((rDc(1)+rD(1,2))/2,(rDc(2)+rD(2,2))/2-0.5,'$\rho_{df}^s$','fontsize',fontSize)
%line joining the first and the last defender
plot([rDc1(1),rDc2(1)],[rDc1(2),rDc2(2)],'b-.')
%x-axis
plot([rDc(1),rDc(1)+8],[rDc(2),rDc(2)],'b--')
phi=atan2(rDc2(2)-rDc1(2),rDc2(1)-rDc1(1))+pi/2;
quiver(rDc(1),rDc(2),cos(phi),sin(phi),4,'color',[0,0,1],'linewidth',1.5,'maxHeadSize',2)
Rp=2;
plot(rDc(1)+Rp*cos(0:phi/100:phi),rDc(2)+Rp*sin(0:phi/100:phi),'b')
Rp=3;
text(rDc(1)+Rp*cos(phi/1.2),rDc(2)+Rp*sin(phi/1.2),'$\phi$','fontsize',fontSize,'color',[0,0,1])

axis(rho_sn*[-1.1,1.1,-1.1,1.1]+[rDc(1),rAcm(1),rDc(2),rAcm(2)])


%angle between the center of defender formation and the ACoM
phi_g=atan2(rAcm(2)-rDc(2),rAcm(1)-rDc(1));
%quiver(rDc(1),rDc(2),cos(phi_g),sin(phi_g),5,'color',[0.5,0,0.5],'maxHeadSize',2)
plot([rDc(1),rAcm(1)],[rDc(2),rAcm(2)],'color',[0.5,0,0.5])
Rp=3.5;
plot(rDc(1)+Rp*cos(0:phi_g/100:phi_g),rDc(2)+Rp*sin(0:phi_g/100:phi_g),'color',[0.5,0,0.5])
Rp=4;
text(rDc(1)+Rp*cos(phi_g/2),rDc(2)+Rp*sin(phi_g/2),'$\theta$','fontsize',fontSize,'color',[0.5,0,0.5])

%plot string lines joining the defenders desired positions
WDString_temp=WDString;
countAPS=0;
for j=1:ND
    for jj=find(WDString_temp(j,:)==1)
        theta=atan2(rD_des(2,jj)-rD_des(2,j),rD_des(1,jj)-rD_des(1,j));
        countAPS=countAPS+1;
        mDD=tan(theta);
        cDD=rD_des(2,j)-mDD*rD_des(1,j);
        rAProjS(1,countAPS)=(mDD*rA(2)+rA(1)-mDD*cDD)/(1+mDD^2);
        rAProjS(2,countAPS)=mDD*rAProjS(1,countAPS)+cDD; 
        thetaProjS(countAPS)=theta+pi/2;
                
        plot([rD_des(1,j)+rhoD*cos(theta),rD_des(1,jj)+rhoD*cos(theta+pi)],[rD_des(2,j)+rhoD*sin(theta),rD_des(2,jj)+rhoD*sin(theta+pi)],'color',[0,0.5,0.5],'linewidth',2)
        plot([rD_des(1,j)+rhoD*cos(theta),rD_des(1,jj)+rhoD*cos(theta+pi)],[rD_des(2,j)+rhoD*sin(theta),rD_des(2,jj)+rhoD*sin(theta+pi)],'w--','linewidth',1)
        %text((rD_des(1,j)+rD_des(1,jj))/2,(rD_des(2,j)+rD_des(2,jj))/2,['$\mathcal{B}_{',num2str(j),'}$'],'fontsize',fontSize)
        %Remove the strings from plot
        WDString_temp(j,jj)=0;
        WDString_temp(jj,j)=0;
    end
end

%text for the desired locations
%text((rD_des(1,1)+rD_des(1,ND))/2,(rD_des(2,1)+rD_des(2,ND))/2,['$\bar{R}_{1}^{',num2str(ND),'} = \bar{R}_d -2b_{d}$'],'fontsize',fontSize)

%Plot radius of the StringNet region
theta=atan2(rD_des(2,ND)-rAcm(2),rD_des(1,ND)-rAcm(1));
R=rho_sn-rhoD;
plot([rAcm(1),rAcm(1)+R*cos(theta)],[rAcm(2),rAcm(2)+R*sin(theta)],'--','color',[0,0.5,0.5])
text([rAcm(1)+rAcm(1)+R*cos(theta)]/2+0,[rAcm(2)+rAcm(2)+R*sin(theta)]/2-1,'$\rho_{sn}$', 'fontsize',fontSize,'color',[0,0.5,0.5])

%line joining the mid-point of first and last defender and center of formation (ACoM)
rDFh=0.5*(rD_des(:,1)+rD_des(:,end));
plot([rAcm(1),rDFh(1,1)],[rAcm(2),rDFh(2,1)],'color',[0,0.5,0.5])
%x-axis
plot([rAcm(1),rAcm(1)+5],[rAcm(2),rAcm(2)],'b--')
phi=atan2(rDFh(2,1)-rAcm(2),rDFh(1,1)-rAcm(1));
%quiver(rAcm(1),rAcm(2),cos(phi),sin(phi),5,'color',[0,0,1],'maxHeadSize',2)
Rp=2;
plot(rAcm(1)+Rp*cos(0:phi/100:phi),rAcm(2)+Rp*sin(0:phi/100:phi),'color',[0,0.5,0.5])
Rp=3;
text(rAcm(1)+Rp*cos(phi/1.2),rAcm(2)+Rp*sin(phi/1.2),'$\phi_{df}^{e*}$','fontsize',fontSize,'color',[0,0.5,0.5])

%Plot the delta agent on the connectivity region
j=5;
RDAcm=norm(rAcm-rD(:,j));
rDProj=(1-rho_Acon/RDAcm)*rAcm+rho_Acon/RDAcm*rD(:,j);

 theta=atan2(rD(2,j)-rAcm(2),rD(2,j)-rAcm(2));
    rTP=[-sin(theta); cos(theta)];
    if rTP'*vD(:,j)<0
        rTP=-rTP;
    end

%Tangent at the projection point
vDProj=rTP'*vD(:,j)*rTP;

%Projection
plot(rDProj(1),rDProj(2),'o','color',colors{1},'markersize',2*markerSize);
text(rDProj(1)-2,rDProj(2)-1.3,'$(\mathbf{r}_{\delta j,c}, \mathbf{v}_{\delta j,c})$','fontsize',fontSize,'color',[0,0,1])

% quiver(rDProj(1),rDProj(2),vDProj(1),vDProj(2),0.1,'linewidth',2.5,'color',colors{1},'MaxHeadSize',.5*headSize)
% text(rDProj(1)-0.3,rDProj(2)-1.8,'$\mathbf{v}_{\delta jk}$','fontsize',fontSize)


%Plot projection of one attacker on the string
j=1;
jj=2;
i=3;
theta=atan2(rD_des(2,jj)-rD_des(2,j),rD_des(1,jj)-rD_des(1,j));
%countAPS=countAPS+1;
mDD=tan(theta);
cDD=rD_des(2,j)-mDD*rD_des(1,j);
rAProjS(1,1)=(mDD*rA(2,i)+rA(1,i)-mDD*cDD)/(1+mDD^2);
rAProjS(2,1)=mDD*rAProjS(1,1)+cDD;

plot(rAProjS(1,1),rAProjS(2,1),'ro','markersize',2*markerSize)
text(rAProjS(1)-3.5,rAProjS(2)+1.3,'$(\mathbf{r}_{\alpha is}, \mathbf{v}_{\alpha is})$','fontsize',fontSize,'color',[1,0,0])


%Annotations
text(0,0,'GatheringPhase','color',[0,0,1], 'fontsize',fontSize)
text(16,-3,'Seeking Phase ','color',[0,0,1], 'fontsize',fontSize)
text(13,-5,'(with superscript s)','color',[0,0,1], 'fontsize',fontSize)
text(21,3.5,'Enclosing Phase','color',[0,0.5,0.5], 'fontsize',fontSize)
%text(-1.5,23,'Herding Phase','color',[0,0.5,0.5], 'fontsize',fontSize)
% text(-3,21.5,'(with center at $\mathbf{r}_{df}^h$)','color',[0,0.5,0.5], 'fontsize',fontSize)

%max footprint that can pass through the 
if(0)
R=rho_sn_max;
theta=pi/2.5;
plot(rAcm(1)+R*circCos,rAcm(2)+R*circSin,'color',colors{4})
plot([rAcm(1),rAcm(1)+R*cos(theta)],[rAcm(2),rAcm(2)+R*sin(theta)],'-','color',colors{4})
de=0.2;
text(de*rAcm(1)+(1-de)*(rAcm(1)+R*cos(theta))-2,de*rAcm(2)+(1-de)*(rAcm(2)+R*sin(theta)),'$\rho_{sn}^{max}$', 'fontsize',fontSize)
R=rho_sn_max-b;
theta=pi/3;
plot(rAcm(1)+R*circCos,rAcm(2)+R*circSin,'--','color',colors{4})
plot([rAcm(1),rAcm(1)+R*cos(theta)],[rAcm(2),rAcm(2)+R*sin(theta)],'--','color',colors{4})
de=0.25;
text(de*rAcm(1)+(1-de)*(rAcm(1)+R*cos(theta))+0.2,de*rAcm(2)+(1-de)*(rAcm(2)+R*sin(theta)),'$\rho_{sn}^{max}-b_d$', 'fontsize',fontSize)
end



%%
%Draw the Two convex curves with their common tangents
colors={[0.5,0.8,0.2],[0,0,1],[0,0,0],[1,0,0.5]};
fontSize=28;
lineWidth=1;
%Plot the obstacles having convex shape
rVO={[0,0;3,0;2,2]',[14,0;15,2;11,2]'};
rho_safe=1.2;
NO=length(rVO);
figure('units','normalized','outerposition',[.1 0.4 .7 .4])
hold all;
%axis equal;
for k=1:NO
    nVO(k)=size(rVO{k},2);
    posVO=rVO{k};
    rCO(:,k)=1/nVO(k)*[sum(posVO,2)];
    posVO=[posVO,posVO(:,1)];
    rVO{k}=posVO;
    %plot(posVO(1,:),posVO(2,:),'color',colors{3})
    
    %Calculate the vertices (outer vertices) of the approximated obstacles
    countV=0;
    for i=1:nVO(k)
        mLi=(posVO(2,i+1)-posVO(2,i))/(posVO(1,i+1)-posVO(1,i));
        if mLi==0
            xi1=posVO(1,i);
            xi2=posVO(1,i);
            yi1=posVO(2,i)+rho_safe;
            yi2=posVO(2,i)-rho_safe;
        elseif isinf(mLi)
            xi1=posVO(1,i)+rho_safe;
            xi2=posVO(1,i)-rho_safe;
            yi1=posVO(2,i);
            yi2=posVO(2,i);
        else
            dx=sqrt(rho_safe^2/(1+1/mLi^2));
            xi1=posVO(1,i)+dx;
            xi2=posVO(1,i)-dx;
            yi1=posVO(2,i)+dx*(-1/mLi);
            yi2=posVO(2,i)-dx*(-1/mLi);
        end
        
        crossProd=cross([posVO(1:2,i+1)-posVO(1:2,i);0],[[xi1;yi1]-posVO(1:2,i);0]);
        countV=countV+1;
        if crossProd(3)<0
            posVO2(:,countV)=[xi1,yi1]';
            if mLi==0
                xip=posVO(1,i+1);
                yip=posVO(2,i+1)+rho_safe;
            elseif isinf(mLi)
                xip=posVO(1,i+1)+rho_safe;
                yip=posVO(2,i+1);
            else
                xip=posVO(1,i+1)+dx;
                yip=posVO(2,i+1)+dx*(-1/mLi);
            end
            
        else
            posVO2(:,countV)=[xi2,yi2]';
            if mLi==0
                xip=posVO(1,i+1);
                yip=posVO(2,i+1)-rho_safe;
            elseif isinf(mLi)
                xip=posVO(1,i+1)-rho_safe;
                yip=posVO(2,i+1);
            else
                xip=posVO(1,i+1)-dx;
                yip=posVO(2,i+1)-dx*(-1/mLi);
            end
        end
        %angVO2(:,countV)=atan2(posVO2(2,countV)-rCO(2,k),posVO2(1,countV)-rCO(1,k));
        
%         if angVO2(:,countV)<0
%             angVO2(:,countV)= angVO2(:,countV)+2*pi;
%         end
%         %corresponding to the next vertex
        countV=countV+1;
        posVO2(:,countV)=[xip,yip]';
        %angVO2(:,countV)=atan2(posVO2(2,countV)-rCO(2,k),posVO2(1,countV)-rCO(1,k));
%         if angVO2(:,countV)<0
%             angVO2(:,countV)= angVO2(:,countV)+2*pi;
%         end
    end
    %Shift the outer vertex corresponding to the first inner vertex to
    %appropriate position
    posVO2=posVO2(:,[countV,1:countV-1]);
    %angVO2=angVO2(:,[countV,1:countV-1]);
    rVO2{k}=posVO2;
    %aVO2{k}=angVO2;
    nV2=countV;
    posVO2=[posVO2,posVO2(:,1)];
    %Find the perimeter and Plot the outer vertices
    PeriO(k)=0;
    for ii=1:nVO(k)
        %Get the inner vertex coordinates
        xV=posVO(1,ii);
        yV=posVO(2,ii);
        AngVOVck(ii,k)=atan2(yV-rCO(2,k),xV-rCO(1,k));
        if  AngVOVck(ii,k)<0 
            AngVOVck(ii,k)= AngVOVck(ii,k)+2*pi;
        end
      
        if ii==1
              AngVOVck1(k)= AngVOVck(ii,k);
        end
        AngVOVck(ii,k)=AngVOVck(ii,k)-AngVOVck1(k);
        if  AngVOVck(ii,k)<0 
            AngVOVck(ii,k)= AngVOVck(ii,k)+2*pi;
        end
        
%         if ii~=1
%             aVO2{k}(2*ii-1)=aVO2{k}(2*ii-1)-aVO2{k}(1);  %shift with respect the first vertex
%             if aVO2{k}(2*ii-1)<0
%                 aVO2{k}(2*ii-1)=aVO2{k}(2*ii-1)+2*pi;
%             end
%             aVO2{k}(2*ii)=aVO2{k}(2*ii)-aVO2{k}(1);  %shift with respect the first vertex
%             if aVO2{k}(2*ii)<0
%                 aVO2{k}(2*ii)=aVO2{k}(2*ii)+2*pi;
%             end
%         else
%             aVO2{k}(2*ii)=aVO2{k}(2*ii)-aVO2{k}(1);  %shift with respect the first vertex
%             if aVO2{k}(2*ii)<0
%                 aVO2{k}(2*ii)=aVO2{k}(2*ii)+2*pi;
%             end
%         end
        
        
        %Find angles corresponding to the outer vertices w.r.t the inner
        %one
        
        ang1=atan2(posVO2(2,2*ii-1)-yV,posVO2(1,2*ii-1)-xV);
        if ang1<0
            ang1=ang1+2*pi;
        end
        
        ang2=atan2(posVO2(2,2*ii)-yV,posVO2(1,2*ii)-xV);
        if ang2<0
            ang2=ang2+2*pi;
        end
        if ang2<ang1
            ang2=ang2+2*pi;
        end
        %Find the perimeters of the segments of the approximated obstacle
        PeriSOk{k}(2*ii-1,1)=rho_safe*(ang2-ang1);
        PeriSOk{k}(2*ii,1)=norm(posVO2(:,2*ii)-posVO2(:,2*ii+1));
        
        AngVO2VOk(2*ii-1,k)=ang1;
        AngVO2VOk(2*ii,k)=ang2;
        
        %plot the circular arc
        plot(xV+rho_safe*cos(ang1:(ang2-ang1)/50:ang2),yV+rho_safe*sin(ang1:(ang2-ang1)/50:ang2),'color',colors{2+k})
        %plot the straight line
        plot(posVO2(1,2*ii:2*ii+1),posVO2(2,2*ii:2*ii+1),'color',colors{2+k});        
    end    
end

ang1=-pi/9;
ang2=3.5*pi/4;
rP1=rVO{1}(:,2)+rho_safe*[cos(ang1),sin(ang1)]';
rP2=rVO{2}(:,3)+rho_safe*[cos(ang2),sin(ang2)]';
rPT1=[cos(ang1+pi/2),sin(ang1+pi/2)]';
rPT2=[cos(ang2+pi/2),sin(ang2+pi/2)]';
plot(rP1(1),rP1(2),'bo')
text(rP1(1)-1.3,rP1(2),'$\gamma_1$','fontsize',fontSize)
plot(rP2(1),rP2(2),'bo')
text(rP2(1)+.3,rP2(2),'$\gamma_2$','fontsize',fontSize)
plot([rP1(1),rP2(1)],[rP1(2),rP2(2)],'color',colors{2})
%Slopes
text(3.5,2.6,'$m_1(\gamma_1)$','fontsize',fontSize,'color',colors{3})
text(9.1,0,'$m_2(\gamma_2)$','fontsize',fontSize,'color',colors{4})
text(5.9,2.1,'$m(\gamma_1,\gamma_2)$','fontsize',fontSize,'color',colors{2})
%x axis
plot([rP1(1),rP1(1)+3],[rP1(2),rP1(2)],'b--')

quiver(rP1(1),rP1(2),rPT1(1),rPT1(2),2.6,'MaxHeadSize',.7*headSize, 'linewidth', lineWidth,'color',colors{3})
quiver(rP2(1),rP2(2),rPT2(1),rPT2(2),2.6,'MaxHeadSize',.7*headSize, 'linewidth', lineWidth,'color',colors{4})



%%
%shortest path 
colors={[0,0,1],[1,0,0],[0,1,0],[0,0,0],[0.4,1,0.1],[0.3,0.3,0.3],[0.5,0.9,0.5],[0,0.5,1],[.1,.51,1],[1,0.2,.61],[0.5,0.5,1],[0.1,0.9,0.1],[0.9,0.4,0.5],[0,0.5,0.4]};
fontSize=18;
global rho_safe;
rho_safe=1;
rVO={[0,1;6,1;6,5;0,5]',[4,-10;11,-10;11,-5;4,-5]'};
fig=figure('units','normalized','outerposition',[.2 0.3 .5 .6])
hold all
axis off
if (0)
tanG=tangentGraph(rVO);
end

lineWidth=1;


%The perimeter of the obstacle
%figure
for k=1:tanG.NO
    posVO2=tanG.rVO2{k};
    posVO2=[posVO2,posVO2(:,1)];
    posVO=tanG.rVO{k};
    rCO(:,k)=1/tanG.nVO(k)*[sum(posVO(:,1:end-1),2)];
    fill(posVO(1,:),posVO(2,:),colors{6})
        text(rCO(1,k)-.5,rCO(2,k),['$\mathcal{O}_',num2str(k),'$'],'color',[1,1,1],'fontsize',fontSize)
    for ii=1:tanG.nVO(k)
        xV=posVO(1,ii);
        yV=posVO(2,ii);
    %Find angles corresponding to the outer vertices w.r.t the inner
        %one
        
        ang1=atan2(posVO2(2,2*ii-1)-yV,posVO2(1,2*ii-1)-xV);
        if ang1<0
            ang1=ang1+2*pi;
        end
        
        ang2=atan2(posVO2(2,2*ii)-yV,posVO2(1,2*ii)-xV);
        if ang2<0
            ang2=ang2+2*pi;
        end
        if ang2<ang1
            ang2=ang2+2*pi;
        end
        %plot the circular arc
        plot(xV+rho_safe*cos(ang1:(ang2-ang1)/50:ang2),yV+rho_safe*sin(ang1:(ang2-ang1)/50:ang2),'color',[0.7,0.7,0.7])
        plot(xV+rho_safe*cos(ang1:(ang2-ang1)/50:ang2),yV+rho_safe*sin(ang1:(ang2-ang1)/50:ang2),'--','color',colors{9})
        plot([xV,xV+rho_safe*cos(ang1)],[yV,yV+rho_safe*sin(ang1)],'--','color',[0.7,0.7,0.7])
        plot([xV,xV+rho_safe*cos(ang2)],[yV,yV+rho_safe*sin(ang2)],'--','color',[0.7,0.7,0.7])
        if k==1
        plot(xV+rho_safe*cos(ang1),yV+rho_safe*sin(ang1),'o','color',colors{9})
         plot(xV+rho_safe*cos(ang2),yV+rho_safe*sin(ang2),'o','color',colors{9})
        end
        %plot the straight line
        plot(posVO2(1,2*ii:2*ii+1),posVO2(2,2*ii:2*ii+1),'color',colors{9}) 
    end
end

%The gamma variables
k=1;
rVO2=tanG.rVO2{k};
text(rVO2(1,1)-3.8,rVO2(2,1),'$\gamma_{\bar{o}1}^1=0$','fontsize',fontSize,'color',colors{9})
text(rVO2(1,2)-0.3,rVO2(2,2)-1,'$\gamma_{\bar{o}1}^2$','fontsize',fontSize,'color',colors{9})
text(rVO2(1,3)-0.3,rVO2(2,3)-1,'$\gamma_{\bar{o}1}^3$','fontsize',fontSize,'color',colors{9})
text(rVO2(1,4)+.3,rVO2(2,4),'$\gamma_{\bar{o}1}^4$','fontsize',fontSize,'color',colors{9})
text(rVO2(1,5)+.3,rVO2(2,5),'$\gamma_{\bar{o}1}^5$','fontsize',fontSize,'color',colors{9})
text(rVO2(1,6)-0.3,rVO2(2,6)+1,'$\gamma_{\bar{o}1}^6$','fontsize',fontSize,'color',colors{9})
text(rVO2(1,7)-0.3,rVO2(2,7)+1,'$\gamma_{\bar{o}1}^7$','fontsize',fontSize,'color',colors{9})
text(rVO2(1,8)-1.8,rVO2(2,8),'$\gamma_{\bar{o}1}^8$','fontsize',fontSize,'color',colors{9})
axis equal

%psi variables
k=1;
rVO=tanG.rVO{k};
kk=1;
ang=atan2(rVO2(2,2*kk-1)-rVO(2,kk),rVO2(1,2*kk-1)-rVO(1,kk));
R=rho_safe/2;
plot([rVO(1,kk),rVO2(1,2*kk-1)],[rVO(2,kk),rVO2(2,2*kk-1)],'--','color',colors{12}')
plot(rVO(1,kk)+R*cos(0:ang/100:ang),rVO(2,kk)+R*sin(0:ang/100:ang),'--','color',colors{12})
R=1.1*R;
text(rVO(1,kk)+.3,rVO(2,kk)+1,['$\psi_{\bar{o}',num2str(k),'}^',num2str(2*kk-1),'$'],'fontsize',fontSize,'color',colors{12})

kk=2;
ang=atan2(rVO2(2,2*kk-1)-rVO(2,kk),rVO2(1,2*kk-1)-rVO(1,kk));
if ang<0
    ang=ang+2*pi;
end
R=rho_safe/2;
plot([rVO(1,kk),rVO2(1,2*kk-1)],[rVO(2,kk),rVO2(2,2*kk-1)],'--','color',colors{12}')
plot(rVO(1,kk)+R*cos(0:ang/100:ang),rVO(2,kk)+R*sin(0:ang/100:ang),'--','color',colors{12})
R=1.1*R;
text(rVO(1,kk)-2,rVO(2,kk)+.8,['$\psi_{\bar{o}',num2str(k),'}^',num2str(2*kk-1),'$'],'fontsize',fontSize,'color',colors{12})

%axes system
quiver(14,5,1,0,3,'color',colors{4},'MaxHeadSize',.7*headSize, 'linewidth', lineWidth)
text(17,4.5,'$\hat{\mathbf{i}}$')
quiver(14,5,0,1,3,'color',colors{4},'MaxHeadSize',.7*headSize, 'linewidth', lineWidth)
text(13.5,8,'$\hat{\mathbf{j}}$')

%Perimeter
k=2;
rVO2=tanG.rVO2{k};
posVO2=rVO2;
posVO=tanG.rVO{k};
plot(rVO2(1,3),rVO2(2,3),'o','color',colors{4})
text(rVO2(1,3)+0.5,rVO2(2,3)-0.5,'$\mathbf{r}_{p1}$','fontsize',fontSize)
plot(rVO2(1,5),rVO2(2,5),'o','color',colors{4})
text(rVO2(1,5)+0.5,rVO2(2,5),'$\mathbf{r}_{p2}$','fontsize',fontSize)
text(12.2,-8,'$Pr(\mathbf{r}_{p1},\mathbf{r}_{p2})$','fontsize',fontSize)
for ii=2
        xV=posVO(1,ii);
        yV=posVO(2,ii);
    %Find angles corresponding to the outer vertices w.r.t the inner one        
        ang1=atan2(posVO2(2,2*ii-1)-yV,posVO2(1,2*ii-1)-xV);
        if ang1<0
            ang1=ang1+2*pi;
        end        
        ang2=atan2(posVO2(2,2*ii)-yV,posVO2(1,2*ii)-xV);
        if ang2<0
            ang2=ang2+2*pi;
        end
        if ang2<ang1
            ang2=ang2+2*pi;
        end
        %plot the circular arc
        plot(xV+rho_safe*cos(ang1:(ang2-ang1)/50:ang2),yV+rho_safe*sin(ang1:(ang2-ang1)/50:ang2),'color',colors{4})        
        %plot the straight line
        plot(posVO2(1,2*ii:2*ii+1),posVO2(2,2*ii:2*ii+1),'color',colors{4}) 
end



%shortest path
rI=[6,-11.5]';
rF=[8,4]';
[Path,tanG_prime,rE_dash] = findShortestPath(tanG,rI,rF);
PlotShortestPath(tanG,Path,tanG_prime,fig.Number,[1,0,0],'-',1)
%Initial and final position text
text(rI(1)+0.5,rI(2),'$\mathbf{r}_0$','color',colors{2},'fontsize',fontSize)
text(rF(1)+0.5,rF(2),'$\mathbf{r}_f$','color',colors{2},'fontsize',fontSize)
%Plot straight line segment joining the points
plot([rI(1),rF(1)],[rI(2),rF(2)],'--','color',colors{10})
text(8,-3,'$L(\mathbf{r}_0,\mathbf{r}_f)$','color',colors{10},'fontsize',fontSize)
%shortest path length
text(-1.0,-3,'$L_s(\mathbf{r}_0,\mathbf{r}_f)$','color',colors{2},'fontsize',fontSize)
%Ellipse filter
Kappa=2;
a=Kappa/2*norm(rI-rF);
c=norm(rI-rF)/2;
b=sqrt(a^2-c^2);
x0=-a:a/50:a;
y0=b/a*sqrt(a^2-x0.^2);
x=[x0,fliplr(x0)];
y=[y0,-fliplr(y0)];
% x1=(rI(1)+rF(1))/2+x;
% y1=(rI(2)+rF(2))/2+y;
thetaIF=atan2(rF(2)-rI(2),rF(1)-rI(1));
Rot=[cos(thetaIF),sin(thetaIF);-sin(thetaIF), cos(thetaIF)];
X=[x;y];
for i=1:length(x)
X2(:,i)=Rot'*X(:,i)+(rI+rF)/2;
end
plot(X2(1,:),X2(2,:),'r--')
text(-4,-20,'$||\mathbf{r}-\mathbf{r}_0||+||\mathbf{r}-\mathbf{r}_f|| =\mathcal{K} ||\mathbf{r}_0-\mathbf{r}_f||$','fontsize',fontSize,'color',colors{2})
quiver(-2,-19.1,1,1,1.8,'color',colors{2},'MaxHeadSize',.7*headSize, 'linewidth', lineWidth)

%velocity labels
