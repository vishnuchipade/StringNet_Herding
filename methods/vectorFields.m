%This script plots different vector fields around circular and rectangular
% ------------------------------
% Author: Vishnu S. Chipade
% Email: vishnuc@umich.edu
% ------------------------------

%obstacles
SetPlotDefaults;
[X,Y]=meshgrid(-20:.5:20,-20:.5:20);
global rS rO rotSense Q G RO RO_max

options = optimset('Display','off','MaxIter',1000);

rI=[0,-10]';
%rS=[-3,39]';
%rO=[6,21;8,-7]';
rS=[6,8]';
rO=[0,0;8,-7]';
NO=size(rO,2);
rotSense=[2,1,2];
RO=2;
%w=sqrt(2)*RO;
w=4;
h=3;
RO_max=6;
for k=1:NO
    xO(k,:)=rO(1,k)+RO*cos(0:pi/100:2*pi);
    yO(k,:)=rO(2,k)+RO*sin(0:pi/100:2*pi);
    ROS(k)=norm(rO(:,k)-rS);
end
P1x=(rS(1)-rI(1))*ones(size(X));
P1y=(rS(2)-rI(2))*ones(size(X));
Q=1;
G=10;
n0=2;
a=w/2*(2^(1/(2*n0)))
b=h/2*(2^(1/(2*n0)))
P2x=zeros(size(X));
P2y=P2x;
P3x=P2x;
P3y=P2x;
P4x=P2x;
P4y=P2x;
for i=1:size(X,1)
    for j=1:size(X,2)
        R2=sqrt((rS(1)-X(i,j))^2+(rS(2)-Y(i,j))^2);
        for k=1:size(rO,2)
            R(k)=sqrt((rO(1,k)-X(i,j))^2+(rO(2,k)-Y(i,j))^2);
            
            RE0=((X(i,j)-rO(1,k))/a)^(2*n0)+((Y(i,j)-rO(2,k))/b)^(2*n0)-1;
            RE(k)=RE0;
            if RE0>0.1
                %                 f1=@(x) x-0.5*(2*(X(i,j)-rO(1,k))/w)^(2/(1-exp(-x)))-0.5*(2*(Y(i,j)-rO(2,k))/h)^(2/(1-exp(-x)))+1;
                %             RE(k)=fsolve(f1,RE0,options);
                %             n=1/(1-exp(-RE(k)));
                n=n0;
                temp=sqrt((b^(2*n)*(abs(X(i,j)-rO(1,1)))^(2*n-1))^2+(a^(2*n)*(abs(Y(i,j)-rO(2,1)))^(2*n-1))^2)
                P4x(i,j)=(-1)^(rotSense(1)+1)*1/RE(k)*a^(2*n)*sign(Y(i,j)-rO(2,1))*(abs(Y(i,j)-rO(2,1)))^(2*n-1)/temp;
                P4y(i,j)=(-1)^(rotSense(1))*1/RE(k)*b^(2*n)*sign(X(i,j)-rO(1,1))*(abs(X(i,j)-rO(1,1)))^(2*n-1)/temp;
            end
            if R(k)>1.1*RO
                flag=1;
            else
                flag=0;
                break;
            end
        end
        for k=1:NO
            if flag==1
                P2x(i,j)=Q*(rS(1)-X(i,j));%/R2^2;
                P2y(i,j)=Q*(rS(2)-Y(i,j));%/R2^2;
                
                P3x(i,j)=P3x(i,j)+(-1)^rotSense(k)*G/(R(k)-RO)*(rO(2,k)-Y(i,j))/R(k)-(-1)^rotSense(k)*G/(ROS(k)-RO)*(rO(2,k)-rS(2))/ROS(k);
                P3y(i,j)=P3y(i,j)+(-1)^(rotSense(k)+1)*G/(R(k)-RO)*(rO(1,k)-X(i,j))/R(k)-(-1)^(rotSense(k)+1)*G/(ROS(k)-RO)*(rO(1,k)-rS(1))/ROS(k);
                
            end
        end
    end
end
%fun=@(x) [Q/((rS(1)-x(1))^2+(rS(2)-x(2))^2)*(rS(1)-x(1))+G/(sqrt((rO(1)-x(1))^2+(rO(2)-x(2))^2)-RO)*(rO(2)-x(2))/sqrt((rO(1)-x(1))^2+(rO(2)-x(2))^2);...
%Q/((rS(1)-x(1))^2+(rS(2)-x(2))^2)*(rS(2)-x(2))-G/(sqrt((rO(1)-x(1))^2+(rO(2)-x(2))^2)-RO)*(rO(1)-x(1))/sqrt((rO(1)-x(1))^2+(rO(2)-x(2))^2)];
fun=@zeroPotential;
[z0]=fsolve(fun,rO(:,1)'+1.1*RO*[1,1])


Px=P1x+P2x+P3x;
Py=P1y+P2y+P3y;

%plot vector fields around circle
if (0)
    figure
    plot(rI(1),rI(2),'*');
    hold on;
    plot(rS(1),rS(2),'rsquare');
    hold on;
    plot(rO(1,1),rO(2,1),'ko')
    hold on;
    for k=1:NO
        plot(xO(k,:),yO(k,:),'k');
        hold on;
    end
    plot(z0(1),z0(2),'ro')
    hold on
    quiver(X,Y,P2x+P3x,P2y+P3y,5)
end

%plot vector field around rectangle
if (0)
    figure
    rectangle('position',[rO(1,1)-w/2,rO(2,1)-h/2,w,h])
    hold on
    %quiver(X,Y,P3x,P3y,2)
    hold on;
    %quiver(X,Y,P4x,P4y,5)
    %quiver(X,Y,P2x+P3x,P2y+P3y,5)
    hold on;
    
    %Super-elliptic vector fields
    xl=10;
    f1=@(x,y) 0.5*(2*(x-rO(1,1))/w)^(2*n)+0.5*(2*(y-rO(2,1))/h)^(2*n)-10;
    ezplot(f1,[-xl,xl])
end
%%
NO=1;
alpha=1;
figure
for k=1:NO
    rectangle('position',[rO(1,k)-w/2,rO(2,k)-h/2,w,h],'FaceColor',[0.5 .5 .5])
    hold on
end
xl=5;
%rectangle('position',[rO(1,2)-w/2,rO(2,2)-h/2,w,h],'FaceColor',[0.5 .5 .5])
colors={'b','r','k','c','m','g'};
for ki=1:5
    K=ki-1+0.3;
    K=K/2;
    %    n=k;
    %   K=log(n/(n-1))
    n=1/(1-exp(-alpha*K));
    %n=3;
    K1(ki)=K;
    N1(ki)=n;
    f2=@(x,y) 0.5*(2*(abs(x-rO(1,1))/w))^(2*n)+0.5*(2*(abs(y-rO(2,1))/h))^(2*n)-1-K;
    ez=ezplot(f2,[-xl,xl]);
    set(ez,'color',colors{ki})
    Legends{ki}=['$E_o$ = ',num2str(K)];
    hold on
    xm=(2*(1+K))^(1/(2*n))*w/2;
    a=w/2*(2^(1/(2*n)));
    b=h/2*(2^(1/(2*n)));
    x0=-xm:xm/9:xm;
    x1=[x0(1),0.50*x0(1)+0.5*x0(2),0.85*x0(1)+0.15*x0(2),x0(2:end-1),0.5*x0(end)+0.5*x0(end-1),0.85*x0(end)+0.15*x0(end-1),x0(end)];
    %x1=xm+(sign(x1).*x1-xm)/(xm^3);
    y1=(2*(1+K)-(2.*(abs(x1./w))).^(2*n)).^(1/(2*n)).*(h/2);
    x1=[x1,x1];
    y1=[y1,-y1];
    for i=1:length(x1)
        temp=sqrt((b^(2*n)*(abs(x1(i)-rO(1,1)))^(2*n-1))^2+(a^(2*n)*(abs(y1(i)-rO(2,1)))^(2*n-1))^2)
        P5x(ki,i)=(-1)^(rotSense(1)+1)*1/K*a^(2*n)*sign(y1(i)-rO(2,1))*(abs(y1(i)-rO(2,1)))^(2*n-1)/temp;
        P5y(ki,i)=(-1)^(rotSense(1))*1/K*b^(2*n)*sign(x1(i)-rO(1,1))*(abs(x1(i)-rO(1,1)))^(2*n-1)/temp;
    end
    
    hold on;
    if (0)
        for k=1:NO
            if ki==1
                quiver(x1+rO(1,k),y1+rO(2,k),P5x(ki,:),P5y(ki,:),.15/K,'color','b');
                hold on;
                % quiver(x1+rO(1,2),y1+rO(2,2),P5x(ki,:),P5y(ki,:),.15/K,'color','b');
            else
                quiver(x1+rO(1,k),y1+rO(2,k),P5x(ki,:),P5y(ki,:),.4/K,'color','b')
                hold on
            end
            %quiver(x1+rO(1,2),y1+rO(2,2),P5x(ki,:),P5y(ki,:),.4/K,'color','b')
        end
    end
    hold on;
end
xlabel('x [m]')
ylabel('y [m]')
title('')
legend(Legends)
% figure
% f2=@(x,y) 0.5*(2*(x-rO(1,1))/w)^(2*n)+0.5*(2*(y-rO(2,1))/h)^(2*n)-1-K;
% ez=ezplot(f2,[-xl,xl])
%%
%super-ellipse around a rectangle
% Kmax=1;
% nmin=1/(1-exp(-alpha*Kmax));
% nmin=2;

clear  X Y
[X,Y]=meshgrid(-10+rO(1,k):.5:10+rO(1,k),-9+rO(2,k):.5:9+rO(2,k));
xl=9;

Ktol=0.05;
Rmax=50;
rho_F=0.581; %formation radius
f1=@(x) x-0.5*(2*(w/2+rho_F)/w)^(2/(1-exp(-x)))-0.5*(2*(h/2+rho_F)/h)^(2/(1-exp(-x)))+1;
K0=fsolve(f1,rho_F+sqrt(w^2+h^2)/2,options);
%K0=0.8;
nmin=1/(1-exp(-alpha*K0));
Ktol=K0;
%Ktol=1.83;
Kmax=10*Ktol
Kmax2=20*Ktol
dKcube=(Kmax2-Kmax)^3;
A=2/dKcube;
B=-3*(Kmax2+Kmax)/dKcube;
C=6*Kmax2*Kmax/dKcube;
D=Kmax2^2*(Kmax2-3*Kmax)/dKcube;
%nmin=1/(1-exp(-0.5));
n=nmin;
a=w/2*(2^(1/(2*n)));
b=h/2*(2^(1/(2*n)));
if(0)
    K=1.1*K0;
    n=1/(1-exp(-alpha*K));
    xm=(2*(1+K))^(1/(2*n))*w/2;
    x1=-xm:xm/9:xm;
    y1=(2*(1+K)-(2.*(abs(x1./w))).^(2*n)).^(1/(2*n)).*(h/2);
    x1=[x1,x1];
    y1=[y1,-y1];
    X=x1;
    Y=y1;
end
Px=zeros(size(X));
Py=Px;
POx=Px;
POy=Py;

betaOS=atan2(rS(2)-rO(2,1),rS(1)-rO(1,k));  %angle between the safe area and the obstacle
beta_barOS=atan2(b^(2*n)*sign(rS(1)-rO(1,1))*(abs(rS(1)-rO(1,1)))^(2*n-1),-a^(2*n)*sign(rS(2)-rO(2,1))*(abs(rS(2)-rO(2,1)))^(2*n-1))  %angle of the tangent at the safe area
% if betaOS<0   % get betaOS between [0,2pi]
%     betaOS=betaOS+2*pi;
% end
% if beta_barOS<0   % get beta_barOS between [0,2pi]
%     beta_barOS=beta_barOS+2*pi;
% end
for i=1:size(X,1)
    for j=1:size(X,2)
        RE0=(abs((X(i,j)-rO(1,k))/a))^(2*nmin)+(abs((Y(i,j)-rO(2,k))/b))^(2*nmin)-1;
        if RE0>=Ktol && RE0<Kmax
            sigma=1;
        elseif RE0>=Kmax && RE0<Kmax2
            sigma=A*RE0^3+B*RE0^2+C*RE0+D;
        else
            sigma=0;
        end
        %RE=RE0;
        betaOF=atan2(Y(i,j)-rO(2,k),X(i,j)-rO(1,k));   %angle between the location and the obstacle
%         if betaOF<0   % get betaOF between [0,2pi]
%             betaOF=betaOF+2*pi;
%         end
        if RE0>=Ktol && RE0<Kmax2
            %                 f1=@(x) x-0.5*(2*(X(i,j)-rO(1,k))/w)^(2/(1-exp(-x)))-0.5*(2*(Y(i,j)-rO(2,k))/h)^(2/(1-exp(-x)))+1;
            %                 RE=fsolve(f1,RE0,options);
            %                 n=1/(1-exp(-RE));
            n=nmin;
            RE=RE0;
            a=w/2*(2^(1/(2*n)));
            b=h/2*(2^(1/(2*n)));
            temp=sqrt((b^(2*n)*(abs(X(i,j)-rO(1,1)))^(2*n-1))^2+(a^(2*n)*(abs(Y(i,j)-rO(2,1)))^(2*n-1))^2);
            beta_bar=atan2(b^(2*n)*sign(X(i,j)-rO(1,1))*(abs(X(i,j)-rO(1,1)))^(2*n-1),-a^(2*n)*sign(Y(i,j)-rO(2,1))*(abs(Y(i,j)-rO(2,1)))^(2*n-1));  %angle of the tangent at the given location
            if beta_bar<0   % get betaOS between [0,2pi]
                beta_bar=beta_bar+2*pi;
            end
             Delta_beta=betaOF-betaOS;
        if Delta_beta<0
            Delta_beta=Delta_beta+2*pi;
        end
        
        Delta_beta_bar=beta_barOS-betaOS;
        if Delta_beta_bar<0
            Delta_beta_bar=Delta_beta_bar+2*pi;
        end
        if Delta_beta<pi
            beta_bar1=beta_bar-Delta_beta_bar+Delta_beta/pi*(Delta_beta_bar-pi);
        else
            beta_bar1=beta_bar-Delta_beta_bar/pi*(Delta_beta-pi);
        end
            Beta1(i,j)=beta_bar1;
            %tangential
            %                  POx(i,j)=(-1)^(rotSense(1)+1)*Rmax/RE*a^(2*n)*sign(Y(i,j)-rO(2,1))*(abs(Y(i,j)-rO(2,1)))^(2*n-1)/temp;
            %             POy(i,j)=(-1)^(rotSense(1))*Rmax/RE*b^(2*n)*sign(X(i,j)-rO(1,1))*(abs(X(i,j)-rO(1,1)))^(2*n-1)/temp;
            %Normalized
            POx(i,j)=(-1)^(rotSense(1))*cos(beta_bar1); %Rmax/RE*
            POy(i,j)=(-1)^(rotSense(1))*sin(beta_bar1); %Rmax/RE*
            mOS=(rS(2)-rO(2,k))/(rS(1)-rO(1,k));
            cOS=mOS*rS(1)+rS(2);
            %symmetric field around the obstacle
            %             if Y(i,j)<=mOS*X(i,j)+cOS
            %                 POx(i,j)=-POx(i,j);
            %                 POy(i,j)=-POy(i,j);
            %             end
            %radial
            %             POy(i,j)=Rmax/RE*a^(2*n)*sign(Y(i,j)-rO(2,k))*(abs(Y(i,j)-rO(2,k)))^(2*n-1)/temp;
            %             POx(i,j)=Rmax/RE*b^(2*n)*sign(X(i,j)-rO(1,k))*(abs(X(i,j)-rO(1,k)))^(2*n-1)/temp;
            
            Px(i,j)=POx(i,j);
            Py(i,j)=POy(i,j);
        end
        if RE0>Ktol
            dist=norm(rS-[X(i,j),Y(i,j)]');
            Px(i,j)=sigma*Px(i,j)+(1-sigma)*(rS(1)-X(i,j))/dist;
            Py(i,j)=sigma*Py(i,j)+(1-sigma)*(rS(2)-Y(i,j))/dist;
        end
    end
end
%%
%f1=@(x) [rS(1)-x(1)+,];
figure

f2=@(x,y) 0.5*(2*(abs(x-rO(1,1)))/w)^(2*nmin)+0.5*(2*(abs(y-rO(2,1)))/h)^(2*nmin)-1-Ktol;
ez=ezplot(f2,[-xl+rO(1,k),xl+rO(1,k),-xl+rO(2,k),xl+rO(2,k)]);
xm=(2*(1+Ktol))^(1/(2*nmin))*w/2;
%text(-xm-1,0,'$E_{o}^m$')
set(ez,'color',[0.5,0,1])
hold on;
f2=@(x,y) 0.5*(2*(abs(x-rO(1,1)))/w)^(2*nmin)+0.5*(2*(abs(y-rO(2,1)))/h)^(2*nmin)-1-Kmax;
ez1=ezplot(f2,[-xl+rO(1,k),xl+rO(1,k),-xl+rO(2,k),xl+rO(2,k)]);
set(ez1,'color',[0.1,1,0])
hold on;
f2=@(x,y) 0.5*(2*(abs(x-rO(1,1)))/w)^(2*nmin)+0.5*(2*(abs(y-rO(2,1)))/h)^(2*nmin)-1-Kmax2;
ez2=ezplot(f2,[-xl+rO(1,k),xl+rO(1,k),-xl+rO(2,k),xl+rO(2,k)]);
set(ez2,'color',[1,.4,0.1])
hold on;
plot(rS(1),rS(2),'ro','markersize',12);
hold on;
legend('$E_o=E_o^m$','$E_o=\bar E_o$','$E_o=E_o^u$','$\textbf{r}_s$','location','northwest')
rectangle('position',[rO(1,1)-w/2,rO(2,1)-h/2,w,h],'FaceColor',[0.5 .5 .5])
hold on;
viscircles([rO(1,1)-w/2,rO(2,1)-h/2;rO(1,1)+w/2,rO(2,1)-h/2;rO(1,1)+w/2,rO(2,1)+h/2;rO(1,1)-w/2,rO(2,1)+h/2],rho_F*ones(1,4),'color','r','LineStyle','--','linewidth',.5)
hold on;
line([rO(1,1)-w/2-rho_F,rO(1,1)+w/2+rho_F],[rO(2,1)-h/2-rho_F,rO(2,1)-h/2-rho_F],'color','r','linestyle','--','linewidth',.5)
hold on
line([rO(2,1)-w/2-rho_F,rO(2,1)-w/2-rho_F],[rO(1,1)-h/2-rho_F,rO(1,1)+h/2+rho_F],'color','r','linestyle','--','linewidth',.5)
hold on;
line([rO(1,1)-w/2-rho_F,rO(1,1)+w/2+rho_F],[rO(2,1)+h/2+rho_F,rO(2,1)+h/2+rho_F],'color','r','linestyle','--','linewidth',.5)
hold on
line([rO(2,1)+w/2+rho_F,rO(2,1)+w/2+rho_F],[rO(1,1)-h/2-rho_F,rO(1,1)+h/2+rho_F],'color','r','linestyle','--','linewidth',.5)
hold on;

quiver(X,Y,Px,Py,.7)
xlabel('x [m]')
ylabel('y [m]')
title('')
%text(rS(1)-1.2,rS(2)+0.95,['$r_s$(',num2str(rS(1)),',',num2str(rS(2)),')'],'fontsize',12,'backgroundcolor',[0.9,0.9,0.9])

if (0)
saveas(gcf,'C:/Data files/Academics@UMich/Research/C-UAS - Vishnu and James/ACC 2019 Latex/figures/vectField1.eps','epsc');
saveas(gcf,'Results/superQuadCont1');
end
