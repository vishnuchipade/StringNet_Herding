% ------------------------------
% Author: Vishnu S. Chipade
% Email: vishnuc@umich.edu
% ------------------------------


%Plots field with multiple obstacles

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
%Ktol=1.83;
pt1=27;
pt2=15;
Kmax=10;
Kmax2=Kmax+30;
dKcube=(Kmax2-Kmax)^3;
A=2/dKcube;
B=-3*(Kmax2+Kmax)/dKcube;
C=6*Kmax2*Kmax/dKcube;
D=Kmax2^2*(Kmax2-3*Kmax)/dKcube;
%nmin=1/(1-exp(-0.5));
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
        sigmaProd=1;
        flagField=1;
        for k=1:size(rO,2)
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
            %Initialize
            Beta_bar1(i,j)=0;
            Beta_bar(i,j)=0;
            if betaOF<0   % get betaOF between [0,2pi]
                betaOF=betaOF+2*pi;
            end
            if RE0>=Ktol && RE0<Kmax2
                %                 f1=@(x) x-0.5*(2*(X(i,j)-rO(1,k))/w)^(2/(1-exp(-x)))-0.5*(2*(Y(i,j)-rO(2,k))/h)^(2/(1-exp(-x)))+1;
                %                 RE=fsolve(f1,RE0,options);
                %                 n=1/(1-exp(-RE));
                %n=nmin;
                RE=RE0;
                %             a=w/2*(2^(1/(2*n)));
                %             b=h/2*(2^(1/(2*n)));
                temp=sqrt((b^(2*n)*(abs(X(i,j)-rO(1,k)))^(2*n-1))^2+(a^(2*n)*(abs(Y(i,j)-rO(2,k)))^(2*n-1))^2);
                beta_bar=atan2(b^(2*n)*sign(X(i,j)-rO(1,k))*(abs(X(i,j)-rO(1,k)))^(2*n-1),-a^(2*n)*sign(Y(i,j)-rO(2,k))*(abs(Y(i,j)-rO(2,k)))^(2*n-1));  %angle of the tangent at the given location
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
                    end
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
                
                Px(i,j)=Px(i,j)+sigma*POx(i,j);
                Py(i,j)=Py(i,j)+sigma*POy(i,j);
                sigmaProd=sigmaProd*(1-sigma);
            end
            if RE0>Ktol
                Px(i,j)=Px(i,j)+sigma*POx(i,j);
                Py(i,j)=Py(i,j)+sigma*POy(i,j);
            else
                Px(i,j)=0;
                Py(i,j)=0;
                flagField=0;
                break;
            end
        end
        if flagField==1
            dist=norm(rS-[X(i,j),Y(i,j)]');
            Px(i,j)=Px(i,j)+sigmaProd*(rS(1)-X(i,j))/dist;
            Py(i,j)=Py(i,j)+sigmaProd*(rS(2)-Y(i,j))/dist;
        end
    end
end
%%
%f1=@(x) [rS(1)-x(1)+,];
if(1)
    figure
    for k=1:size(rO,2)
        xsup=a*(Ktol+1)^(1/(2*n))*cos(0:pi/50:pi/2).^(1/n);
        ysup=b*(Ktol+1)^(1/(2*n))*sin(0:pi/50:pi/2).^(1/n);
        xsup=[xsup,-flip(xsup),-xsup,flip(xsup)];
        ysup=[ysup,flip(ysup),-ysup,-flip(ysup)];
        xsup=xsup+rO(1,k);
        ysup=ysup+rO(2,k);
        plot(xsup,ysup,'color',[0.5,0,1],'linewidth',0.5)
        %set(ez,'color',[0.5,0,1])
        hold on;
        xsup=a*(Kmax+1)^(1/(2*n))*cos(0:pi/50:pi/2).^(1/n);
        ysup=b*(Kmax+1)^(1/(2*n))*sin(0:pi/50:pi/2).^(1/n);
        xsup=[xsup,-flip(xsup),-xsup,flip(xsup)];
        ysup=[ysup,flip(ysup),-ysup,-flip(ysup)];
        xsup=xsup+rO(1,k);
        ysup=ysup+rO(2,k);
        plot(xsup,ysup,'color',[0.3,.91,0],'linewidth',0.5)
        hold on;
        xsup=a*(Kmax2+1)^(1/(2*n))*cos(0:pi/50:pi/2).^(1/n);
        ysup=b*(Kmax2+1)^(1/(2*n))*sin(0:pi/50:pi/2).^(1/n);
        xsup=[xsup,-flip(xsup),-xsup,flip(xsup)];
        ysup=[ysup,flip(ysup),-ysup,-flip(ysup)];
        xsup=xsup+rO(1,k);
        ysup=ysup+rO(2,k);
        plot(xsup,ysup,'color',[0,.0,0],'linewidth',0.5)
        hold on;
        plot(rS(1),rS(2),'ro','markersize',12);
        hold on;
        %legend(['$E_o=E_o^m=$',num2str(Ktol,3)],['$E_o=\bar E_o=$',num2str(Kmax,3)],['$E_o=E_o^u=$',num2str(Kmax2,3)],'$\textbf{r}_s$','location','northwest')
        rectangle('position',[rO(1,k)-w/2,rO(2,k)-h/2,w,h],'FaceColor',[0.5 .5 .5])
        hold on;
        viscircles([rO(1,k)-w/2,rO(2,k)-h/2;rO(1,k)+w/2,rO(2,k)-h/2;rO(1,k)+w/2,rO(2,k)+h/2;rO(1,k)-w/2,rO(2,k)+h/2],rho_F*ones(1,4),'color','r','LineStyle','--','linewidth',.5)
        hold on;
        line([rO(1,k)-w/2-rho_F,rO(1,k)+w/2+rho_F],[rO(2,k)-h/2-rho_F,rO(2,k)-h/2-rho_F],'color','r','linestyle','--','linewidth',.5)
        hold on
        line([rO(1,k)-w/2-rho_F,rO(1,k)-w/2-rho_F],[rO(2,k)-h/2-rho_F,rO(2,k)+h/2+rho_F],'color','r','linestyle','--','linewidth',.5)
        hold on;
        line([rO(1,k)-w/2-rho_F,rO(1,k)+w/2+rho_F],[rO(2,k)+h/2+rho_F,rO(2,k)+h/2+rho_F],'color','r','linestyle','--','linewidth',.5)
        hold on
        line([rO(1,k)+w/2+rho_F,rO(1,k)+w/2+rho_F],[rO(2,k)-h/2-rho_F,rO(2,k)+h/2+rho_F],'color','r','linestyle','--','linewidth',.5)
        hold on;
        hold on;
    end
    
    hf=quiver(X,Y,Px,Py,.7)
    set(hf,'color',[.1,.1,1.0])
    xlabel('x [m]')
    ylabel('y [m]')
    title('')
    %text(rS(1)-1.2,rS(2)+0.95,['$r_s$(',num2str(rS(1)),',',num2str(rS(2)),')'],'fontsize',12,'backgroundcolor',[0.9,0.9,0.9])
    
end

% headSize=1.1;
% fontSize=16;
% arrowLength=3;
% plot(rS(1),rS(2),'ro','markersize',8);
% hold on;
% text(rS(1)+0.52,rS(2)+0.15,['$\textbf{r}_s$(',num2str(rS(1)),', ',num2str(rS(2)),')'],'fontsize',12)
% hold on;
% %Tangent and field line at safe location
% quiver(rS(1),rS(2),cos(betaOS)/1.2,sin(betaOS)/1.2,arrowLength,'linewidth',2,'color',[1 0 0],'MaxHeadSize',0.9*headSize)
% quiver(rS(1),rS(2),cos(beta_barOS),sin(beta_barOS),arrowLength,'linewidth',2,'color',[1 0 .7],'MaxHeadSize',0.9*headSize)
% rd=0.8;
% dbeta_barOS=beta_barOS-betaOS;
% plot(rS(1)+rd*cos(betaOS:dbeta_barOS/10:beta_barOS),rS(2)+rd*sin(betaOS:dbeta_barOS/10:beta_barOS),'color',[1 0 0.7])
% rd=2.3;
% text(rS(1)+rd*cos(beta_barOS/1.2),rS(2)+rd*sin(beta_barOS/1.2),'\boldmath$\Delta\bar\beta_{ok}$','fontsize',fontSize,'color',[1 0 0.7])
% %the obstacle
% rectangle('position',[rO(1,1)-w/2,rO(2,1)-h/2,w,h])%,'FaceColor',[0.5 .5 .5])
% hold on;
% text(rO(1,1)-0.15,rO(2,1)-0.4,['$\textbf{r}_{ok}$(',num2str(rO(1,1)),', ',num2str(rO(2,1)),')'],'fontsize',12)
% text(rO(1,1)-1.9,rO(2,1)+1,'$\mathcal{O}_{k}$','fontsize',fontSize)
% hold on;
% f2=@(x,y) 0.5*(2*(abs(x-rO(1,1)))/w)^(2*nmin)+0.5*(2*(abs(y-rO(2,1)))/h)^(2*nmin)-1-RE;
% ez=ezplot(f2,[-xm+rO(1,k),xm+rO(1,k),-xm+rO(2,k),xm+rO(2,k)]);
% set(ez,'linewidth',2,'color',[0.1,0.8,0.6])
% hold on;
% 
% 
% text(X(pt1,pt2)-2.02,Y(pt1,pt2)-.80,'$\textbf{r}(x,y)$','fontsize',12);
% hold on;
% %tangent
% line([X(pt1,pt2),3+X(pt1,pt2)],[Y(pt1,pt2),Y(pt1,pt2)],'linestyle','--','color','b')
% hold on;
% quiver(X(pt1,pt2),Y(pt1,pt2),cos(Beta_bar(pt1,pt2)),sin(Beta_bar(pt1,pt2)),arrowLength,'linewidth',2,'color',[0.1 0.7 0.1],'MaxHeadSize',headSize)
% hold on
% rd=1.2;
% plot(X(pt1,pt2)+rd*cos(0:Beta_bar(pt1,pt2)/10:Beta_bar(pt1,pt2)),Y(pt1,pt2)+rd*sin(0:Beta_bar(pt1,pt2)/10:Beta_bar(pt1,pt2)),'color',[0.1 0.7 0.1])
% hold on;
% rd=2.7;
% text(X(pt1,pt2)+rd*cos(Beta_bar(pt1,pt2)/1.30),Y(pt1,pt2)+rd*sin(Beta_bar(pt1,pt2)/1.30),'\boldmath$\bar{{\beta}}_{ok}(\textbf{r})$','fontsize',fontSize,'color',[0.1 0.7 0.1])
% %vector field
% quiver(X(pt1,pt2),Y(pt1,pt2),Px(pt1,pt2),Py(pt1,pt2),arrowLength,'linewidth',2,'color',[1 0.5 0],'MaxHeadSize',headSize)
% hold on;
% rd=1.8;
% plot(X(pt1,pt2)+rd*cos(0:Beta_bar1(pt1,pt2)/10:Beta_bar1(pt1,pt2)),Y(pt1,pt2)+rd*sin(0:Beta_bar1(pt1,pt2)/10:Beta_bar1(pt1,pt2)),'color',[1 0.5 0])
% hold on
% rd=2.2;
% text(X(pt1,pt2)+rd*cos(Beta_bar1(pt1,pt2)/2.6),Y(pt1,pt2)+rd*sin(Beta_bar1(pt1,pt2)/2.6),'\boldmath$\bar\beta_{ok}^{\prime}(\textbf{r})$','fontsize',fontSize,'color',[1 0.5 0])
% %current location
% line([rO(1,1),7+rO(1,1)],[rO(2,1),rO(2,1)],'linestyle','--','color','b')
% hold on
% plot([rO(1,1),X(pt1,pt2)],[rO(2,1),Y(pt1,pt2)],'b')
% hold on;
% rd=1.15;
% plot(rO(1,1)+rd*cos(0:BetaOF(pt1,pt2)/10:BetaOF(pt1,pt2)),rO(2,1)+rd*sin(0:BetaOF(pt1,pt2)/10:BetaOF(pt1,pt2)),'b')
% hold on;
% rd=1.5;
% text(rO(1,1)+rd*cos(BetaOF(pt1,pt2)/5),rO(2,1)+rd*sin(BetaOF(pt1,pt2)/5),'\boldmath$\beta_{ok}(\textbf{r})$','fontsize',fontSize,'color','b')
% %safe location
% l=line([rO(1,1),rS(1)],[rO(2,1),rS(2)],'linestyle','--','color','r')
% hold on
% rd=4.5;
% plot(rO(1,1)+rd*cos(0:betaOS/10:betaOS),rO(2,1)+rd*sin(0:betaOS/10:betaOS),'r')
% hold on;
% rd=4.8;
% text(rO(1,1)+rd*cos(betaOS/2),rO(2,1)+rd*sin(betaOS/2),'\boldmath$\beta_{ok}(\textbf{r}_s)$','fontsize',fontSize,'color','r')
% hold on
% rd=2.65;
% plot(rO(1,1)+rd*cos(BetaOF(pt1,pt2):(betaOS-BetaOF(pt1,pt2))/10:betaOS),rO(2,1)+rd*sin(BetaOF(pt1,pt2):(betaOS-BetaOF(pt1,pt2))/10:betaOS),'r')
% hold on;
% rd=3.3;
% text(rO(1,1)+rd*cos(BetaOF(pt1,pt2)/1.1),rO(2,1)+rd*sin(BetaOF(pt1,pt2)/1.1),'\boldmath$\Delta\beta_{ok}$','fontsize',fontSize,'color','r')
% hold on
% 
% xlabel('x [m]')
% ylabel('y [m]')
% title('')
% xlim([-1*xl,max(xl,1.1*rS(1))])
% ylim([-0.6*xl,1.3*rS(2)])
