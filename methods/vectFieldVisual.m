% ------------------------------
% Author: Vishnu S. Chipade
% Email: vishnuc@umich.edu
% ------------------------------


options = optimset('Display','off','MaxIter',1000);

rI=[0,-10]';
%rS=[-3,39]';
%rO=[6,21;8,-7]';
rS=[2,6]';
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
% for k=1:NO
%     xO(k,:)=rO(1,k)+RO*cos(0:pi/100:2*pi);
%     yO(k,:)=rO(2,k)+RO*sin(0:pi/100:2*pi);
%     ROS(k)=norm(rO(:,k)-rS);
% end
[X,Y]=meshgrid(-10:.5:10,-9:.5:9);
xl=19;
alpha=1;



Ktol=0.2;
%Ktol=1.83;
Kmax=10*Ktol
%K=4.3/2;
Kmax2=20*Ktol
dKcube=(Kmax2-Kmax)^3;
A=2/dKcube;
B=-3*(Kmax2+Kmax)/dKcube;
C=6*Kmax2*Kmax/dKcube;
D=Kmax2^2*(Kmax2-3*Kmax)/dKcube;


%Plot the generic super-elliptic contours
Karr=[0.25,.75,1.75,2.25]
figure
colors={'b',[0.6,0.1,.3],'k',[0.1,0.8,0.6],'m','g'};
for ki=1:4
    K=ki-1+0.5;
    K=Karr(ki)
    %    n=k;
    %   K=log(n/(n-1))
    n=1/(1-exp(-alpha*K));
    %n=3;
    K1(ki)=K;
    N1(ki)=n;
    f2=@(x,y) 0.5*(2*(abs(x-rO(1,1))/w))^(2*n)+0.5*(2*(abs(y-rO(2,1))/h))^(2*n)-1-K;
    ez=ezplot(f2,[-xl,xl]);
    set(ez,'color',colors{ki},'LineStyle','-.')
    Legends{ki}=['$E_o$ = ',num2str(K)];
    hold on
    xm=(2*(1+K))^(1/(2*n))*w/2;
    a=w/2*(2^(1/(2*n)));
    b=h/2*(2^(1/(2*n)));
    x0=-xm:xm/9:xm;
    x1=[x0(1),0.50*x0(1)+0.5*x0(2),0.85*x0(1)+0.15*x0(2),x0(2:end-1),0.5*x0(end)+0.5*x0(end-1),0.85*x0(end)+0.15*x0(end-1),x0(end)];
    %x1=xm+(sign(x1).*x1-xm)/(xm^3);
    y1=(abs(2*(1+K)-(2.*(abs(x1./w))).^(2*n))).^(1/(2*n)).*(h/2);
    x1=[x1,x1];
    y1=[y1,-y1];
    for i=1:length(x1)
        temp=sqrt((b^(2*n)*(abs(x1(i)-rO(1,1)))^(2*n-1))^2+(a^(2*n)*(abs(y1(i)-rO(2,1)))^(2*n-1))^2);
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
legend(Legends,'Location','northwest')



%%Plot the main contour with angles defined
if(1)
    n=1/(1-exp(-alpha*K));
    nmin=n;
    a=w/2*(2^(1/(2*n)));
    b=h/2*(2^(1/(2*n)));
    xm=(2*(1+K))^(1/(2*n))*w/2;
    x0=-xm:xm/9:xm;
     x1=[x0(1),0.50*x0(1)+0.5*x0(2),0.85*x0(1)+0.15*x0(2),x0(2:end-1),0.5*x0(end)+0.5*x0(end-1),0.85*x0(end)+0.15*x0(end-1),x0(end)];
    for i=1:length(x1)
        y1(i)=(abs(2*(1+K)-(2*(abs(x1(i)/w)))^(2*n)))^(1/(2*n))*(h/2);
    end
   
    x1=[x1,x1];
    y1=[y1,-y1];
    X=x1;
    Y=y1;
end
%X=x1;
%Y=y1;
Px=zeros(size(X));
Py=Px;
POx=Px;
POy=Py;

betaOS=atan2(rS(2)-rO(2,1),rS(1)-rO(1,1));  %angle between the safe area and the obstacle
beta_barOS=atan2(b^(2*n)*sign(rS(1)-rO(1,1))*(abs(rS(1)-rO(1,1)))^(2*n-1),-a^(2*n)*sign(rS(2)-rO(2,1))*(abs(rS(2)-rO(2,1)))^(2*n-1))  %angle of the tangent at the safe area
if betaOS<0   % get betaOS between [0,2pi]
    betaOS=betaOS+2*pi;
end
if beta_barOS<0   % get beta_barOS between [0,2pi]
    beta_barOS=beta_barOS+2*pi;
end
for i=1:size(X,1)
    for j=1:size(X,2)
        RE0=(abs((X(i,j)-rO(1,k))/a))^(2*n)+(abs((Y(i,j)-rO(2,k))/b))^(2*n)-1;
        if RE0>=Ktol && RE0<Kmax
            sigma=1;
        elseif RE0>=Kmax && RE0<Kmax2
            sigma=A*RE0^3+B*RE0^2+C*RE0+D;
        else
            sigma=0;
        end
        %RE=RE0;
        betaOF=atan2(Y(i,j)-rO(2,k),X(i,j)-rO(1,k));   %angle between the location and the obstacle
        BetaOF(i,j)=betaOF;
        if betaOF<0   % get betaOF between [0,2pi]
            betaOF=betaOF+2*pi;
        end
        if RE0>=Ktol && RE0<Kmax2
            %                 f1=@(x) x-0.5*(2*(X(i,j)-rO(1,k))/w)^(2/(1-exp(-x)))-0.5*(2*(Y(i,j)-rO(2,k))/h)^(2/(1-exp(-x)))+1;
            %                 RE=fsolve(f1,RE0,options);
            %                 n=1/(1-exp(-RE));
            %n=nmin;
            RE=RE0;
            a=w/2*(2^(1/(2*n)));
            b=h/2*(2^(1/(2*n)));
            temp=sqrt((b^(2*n)*(abs(X(i,j)-rO(1,1)))^(2*n-1))^2+(a^(2*n)*(abs(Y(i,j)-rO(2,1)))^(2*n-1))^2);
            beta_bar=atan2(b^(2*n)*sign(X(i,j)-rO(1,1))*(abs(X(i,j)-rO(1,1)))^(2*n-1),-a^(2*n)*sign(Y(i,j)-rO(2,1))*(abs(Y(i,j)-rO(2,1)))^(2*n-1));  %angle of the tangent at the given location
            if beta_bar<0   % get betaOS between [0,2pi]
                beta_bar=beta_bar+2*pi;
            end
            %beta=atan2(-b^(2*n)*(abs(X(i,j)-rO(1,1)))^(2*n-1),a^(2*n)*(abs(Y(i,j)-rO(2,1)))^(2*n-1));
            if betaOS<pi
                if betaOF<=betaOS
                    beta_bar1=beta_bar-((-betaOS+pi+betaOF)/pi)*(beta_barOS-betaOS);
                elseif betaOF>betaOS+pi
                    % beta_bar1=betaOS+((betaOS-beta0+pi)/pi)*(beta0-betaOS);
                    beta_bar1=beta_bar+((betaOS+pi-betaOF)/pi)*(beta_barOS-betaOS);
                else
                    % beta_bar1=betaOS+((beta0-betaOS)/pi)*(betaOF-betaOS);
                    
                    beta_bar1=beta_bar-(beta_barOS-betaOS)+(beta_barOS-betaOS-pi)*(betaOF-betaOS)/pi;
                end
            elseif betaOS>=pi && betaOS<1.5*pi
                if betaOF<=betaOS && betaOF>betaOS-pi
                    beta_bar1=beta_bar+((betaOS-betaOF-pi)/pi)*(beta_barOS-betaOS);
                elseif betaOF>betaOS
                    beta_bar1=beta_bar-(beta_barOS-betaOS)+(beta_barOS-betaOS-pi)*(betaOF-betaOS)/pi;
                else
                    beta_bar1=beta_bar-(beta_barOS-betaOS)+(beta_barOS-betaOS-pi)*(betaOF+2*pi-betaOS)/pi;
                end
            else
                if beta_barOS<pi/2
                beta_barOS=beta_barOS+2*pi;
                end;
                if betaOF<=betaOS && betaOF>betaOS-pi
                    beta_bar1=beta_bar+((betaOS-betaOF-pi)/pi)*(beta_barOS-betaOS);
                elseif betaOF>betaOS
                    beta_bar1=beta_bar-(beta_barOS-betaOS)+(beta_barOS-betaOS-pi)*(betaOF-betaOS)/pi;
                else
                    beta_bar1=beta_bar-(beta_barOS-betaOS)+(beta_barOS-betaOS-pi)*(betaOF+2*pi-betaOS)/pi;
                end
            end
            Beta_bar1(i,j)=beta_bar1;
            Beta_bar(i,j)=beta_bar;
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
%figure
%the safe location
headSize=1;
fontSize=15;
plot(rS(1),rS(2),'ro','markersize',8);
hold on;
text(rS(1)+0.22,rS(2)+0.05,['$\textbf{r}_s$(',num2str(rS(1)),', ',num2str(rS(2)),')'],'fontsize',12)
hold on;
quiver(rS(1),rS(2),cos(betaOS)/1.2,sin(betaOS)/1.2,2,'linewidth',2,'color',[1 0 0],'MaxHeadSize',0.9*headSize)
quiver(rS(1),rS(2),cos(beta_barOS),sin(beta_barOS),2,'linewidth',2,'color',[1 0 .7],'MaxHeadSize',0.9*headSize)
rd=0.7;
dbeta_barOS=beta_barOS-betaOS;
plot(rS(1)+rd*cos(betaOS:dbeta_barOS/10:beta_barOS),rS(2)+rd*sin(betaOS:dbeta_barOS/10:beta_barOS),'color',[1 0 0.7])
rd=1.7;
text(rS(1)+rd*cos(beta_barOS/1.2),rS(2)+rd*sin(beta_barOS/1.2),'$\Delta\bar\beta_{ok}$','fontsize',fontSize,'color',[1 0 0.7])
%the obstacle
rectangle('position',[rO(1,1)-w/2,rO(2,1)-h/2,w,h])%,'FaceColor',[0.5 .5 .5])
hold on;
text(rO(1,1)-0.15,rO(2,1)-0.4,['$\textbf{r}_{ok}$(',num2str(rO(1,1)),', ',num2str(rO(2,1)),')'],'fontsize',12)
text(rO(1,1)-1.9,rO(2,1)+1,'$\mathcal{O}_{k}$','fontsize',fontSize)
hold on;
f2=@(x,y) 0.5*(2*(abs(x-rO(1,1)))/w)^(2*nmin)+0.5*(2*(abs(y-rO(2,1)))/h)^(2*nmin)-1-RE;
ez=ezplot(f2,[-xm+rO(1,k),xm+rO(1,k),-xm+rO(2,k),xm+rO(2,k)]);
set(ez,'linewidth',2,'color',[0.1,0.8,0.6])
hold on;
pt=length(X)/2-5

text(X(pt)+0.22,Y(pt)-.32,'$\textbf{r}(x,y)$','fontsize',12);
hold on;
%tangent
line([X(pt),5],[Y(pt),Y(pt)],'linestyle','--','color','b')
hold on;
quiver(X(pt),Y(pt),cos(Beta_bar(pt)),sin(Beta_bar(pt)),2,'linewidth',2,'color',[0.1 0.7 0.1],'MaxHeadSize',headSize)
hold on
rd=0.5;
plot(X(pt)+rd*cos(0:Beta_bar(pt)/10:Beta_bar(pt)),Y(pt)+rd*sin(0:Beta_bar(pt)/10:Beta_bar(pt)),'color',[0.1 0.7 0.1])
hold on;
rd=1.5;
text(X(pt)+rd*cos(Beta_bar(pt)/1.15),Y(pt)+rd*sin(Beta_bar(pt)/1.15),'$\bar\beta_{ok}(\textbf{r})$','fontsize',fontSize,'color',[0.1 0.7 0.1])
%vector field
quiver(X(pt),Y(pt),Px(pt),Py(pt),2,'linewidth',2,'color',[1 0.5 0],'MaxHeadSize',headSize)
hold on;
rd=1;
plot(X(pt)+rd*cos(0:Beta_bar1(pt)/10:Beta_bar1(pt)),Y(pt)+rd*sin(0:Beta_bar1(pt)/10:Beta_bar1(pt)),'color',[1 0.5 0])
hold on
rd=1.2;
text(X(pt)+rd*cos(Beta_bar1(pt)/1.63),Y(pt)+rd*sin(Beta_bar1(pt)/1.63),'$\bar\beta_{ok}^{\prime}(\textbf{r})$','fontsize',fontSize,'color',[1 0.5 0])
%current location
line([rO(1,1),5+rO(1,1)],[rO(2,1),rO(2,1)],'linestyle','--','color','b')
hold on
plot([rO(1,1),X(pt)],[rO(2,1),Y(pt)],'b')
hold on;
rd=2.65;
plot(rO(1,1)+rd*cos(0:BetaOF(pt)/10:BetaOF(pt)),rO(2,1)+rd*sin(0:BetaOF(pt)/10:BetaOF(pt)),'b')
hold on;
rd=2.89;
text(rO(1,1)+rd*cos(BetaOF(pt)/1.73),rO(2,1)+rd*sin(BetaOF(pt)/1.73),'$\beta_{ok}(\textbf{r})$','fontsize',fontSize,'color','b')
%safe location
l=line([rO(1,1),rS(1)],[rO(2,1),rS(2)],'linestyle','--','color','r')
hold on
rd=.8;
plot(rO(1,1)+rd*cos(0:betaOS/10:betaOS),rO(2,1)+rd*sin(0:betaOS/10:betaOS),'r')
hold on;
rd=1;
text(rO(1,1)+rd*cos(betaOS/2.5),rO(2,1)+rd*sin(betaOS/2.5),'$\beta_{ok}(\textbf{r}_s)$','fontsize',fontSize,'color','r')
hold on
rd=2.0;
plot(rO(1,1)+rd*cos(BetaOF(pt):(betaOS-BetaOF(pt))/10:betaOS),rO(2,1)+rd*sin(BetaOF(pt):(betaOS-BetaOF(pt))/10:betaOS),'r')
hold on;
rd=2.6
text(rO(1,1)+rd*cos(betaOS/1.05),rO(2,1)+rd*sin(betaOS/1.05),'$\Delta\beta_{ok}$','fontsize',fontSize,'color','r')
hold on

xlabel('x [m]')
ylabel('y [m]')
title('')
xlim([-1.1*xm,max(1.2*xm,1.5*rS(1))])
ylim([-0.9*xm,1.3*rS(2)])






if (0)
saveas(gcf,'C:/Data files/Academics@UMich/Research/C-UAS - Vishnu and James/ACC 2019 Latex/figures/vectField1.eps','epsc');
saveas(gcf,'Results/superQuadCont1');
end
