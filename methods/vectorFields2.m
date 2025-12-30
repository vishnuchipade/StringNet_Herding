% ------------------------------
% Author: Vishnu S. Chipade
% Email: vishnuc@umich.edu
% ------------------------------

options = optimset('Display','off','MaxIter',1000);

AllParameters

%rS=rP;
ROS=norm(rS-rO(1:2,1));
NO=size(rO,2);
rotSense=[2,1,2];
RO_max=6;
k=1;

clear  X Y
xl=150;
[X,Y]=meshgrid(-xl+rO(1,1):6:xl+rO(1,1),-xl+rO(2,1):6:xl+rO(2,1));

alpha=1;

n=nO(k);
a=aO(k);
b=bO(k); 

rho_F=0.5;

pt1=27;
pt2=15;

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
ii=1;
for i=1:size(X,1)
    for j=1:size(X,2)
        rA=[X(i,j);Y(i,j)];
        rA_goal=rS;
        vA=[0;0];
        vA_goal=vA;
        
        F_A=[0,0]';
        F_A_dot=[0,0]';
            R_OA=norm(rA-rO(1:2,k));
            R_AgoalA=norm(rA_goal-rA);
             EmAO=E_m_AO(ii,k);
            EbarAO=E_bar_AO(ii,k);
            EuAO=E_u_AO(ii,k);
            EvAO=E_u_AO(ii,k);
            E_AOk=(abs((rA(1)-rO(1,k))/a))^(2*n)+(abs((rA(2)-rO(2,k))/b))^(2*n)-1;
            E_AOk_dot=2*n*[(abs((rA(1)-rO(1,k))/a))^(2*n)/(rA(1)-rO(1,k)), (abs((rA(2)-rO(2,k))/b))^(2*n)/(rA(2)-rO(2,k))]*vA;
            %Modify the boundary elliptic distances
            E_Ades_Ok=(abs((rA_goal(1)-rO(1,k))/a))^(2*n)+(abs((rA_goal(2)-rO(2,k))/b))^(2*n)-1;
            if E_Ades_Ok<EuAO
                EuAO=E_Ades_Ok;
                EbarAO=(E_Ades_Ok-EmAO)/2;
            end
            arr_EAO(ii,k)=E_AOk;
            if E_AOk>EmAO &&  E_AOk<EbarAO
                sigma=1;   
                sigma_dot=0;
                sigma_bar=1;
            elseif  E_AOk>EbarAO &&  E_AOk<EuAO
                sigma=A_A_O(ii,k)*E_AOk^3+B_A_O(ii,k)*E_AOk^2+C_A_O(ii,k)*E_AOk+D_A_O(ii,k);
                sigma_dot=(3*A_A_O(ii,k)*E_AOk^2+2*B_A_O(ii,k)*E_AOk+C_A_O(ii,k))*E_AOk_dot;
                sigma_bar=1;
            elseif E_AOk>EuAO &&  E_AOk<EvAO
                sigma_bar=A_bar_A_O(ii,k)*E_AOk^3+B_bar_A_O(ii,k)*E_AOk^2+C_bar_A_O(ii,k)*E_AOk+D_bar_A_O(ii,k);
                sigma=0;   
                sigma_dot=0;
            else
                sigma=0;   
                sigma_dot=0;
                sigma_bar=0;
            end
            Sigma(k)=sigma;
            Sigma_dot(k)=sigma_dot;
            
            %find vector field direction
            betaOS=atan2(rA_goal(2)-rO(2,k),rA_goal(1)-rO(1,k));  %angle between the desired location and the obstacle
            betaOA=atan2(rA(2)-rO(2,k),rA(1)-rO(1,k));   %angle between the location and the obstacle
            betaOA_dot=[-(rA(2)-rO(2,k)), (rA(1)-rO(1,k))]*vA/R_OA^2;
            beta_barOS=atan2(b^(2*n)*sign(rA_goal(1)-rO(1,k))*(abs(rA_goal(1)-rO(1,k)))^(2*n-1),-a^(2*n)*sign(rA_goal(2)-rO(2,k))*(abs(rA_goal(2)-rO(2,k)))^(2*n-1));  %angle of the tangent at the safe area
            betaSA=atan2(rA(2)-rA_goal(2),rA(1)-rA_goal(1));
            betaSA_dot=[-(rA(2)-rA_goal(2)), (rA(1)-rA_goal(1))]*(vA-vA_goal)/R_AgoalA^2;
            beta_bar=atan2(b^(2*n)*sign(rA(1)-rO(1,k))*(abs(rA(1)-rO(1,k)))^(2*n-1),-a^(2*n)*sign(rA(2)-rO(2,k))*(abs(rA(2)-rO(2,k)))^(2*n-1));  %angle of the tangent at the given location
            if beta_bar<0   % get betaOS between [0,2pi]
                beta_bar=beta_bar+2*pi;
            end
            Delta_beta=betaOA-betaOS;
            betaA0=pi+pi/4;
            Delta_beta=betaOA-betaA0;
            if Delta_beta<0
                Delta_beta=Delta_beta+2*pi;
            end
            
            Delta_beta_bar=beta_barOS-betaOS;
            if Delta_beta_bar<0
                Delta_beta_bar=Delta_beta_bar+2*pi;
            end
            
             crossProd=cross([cos(betaA0);sin(betaA0);0],[cos(beta_bar);sin(beta_bar);0]);
            if crossProd(3)>0
                beta_bar1=betaA0;
                beta_bar1_dot=betaSA_dot;            
            else
                if Delta_beta>pi
%                 beta_bar1=beta_bar-Delta_beta_bar+Delta_beta/pi*(Delta_beta_bar-pi);
%                 beta_bar1_dot=(sin(2*betaOA)/(0.5*((cos(2*betaOA))^2+1)) + (Delta_beta_bar-pi)/pi)*betaOA_dot;
                beta_bar1=beta_bar;
                 else
%                 beta_bar1=beta_bar-Delta_beta_bar/pi*(Delta_beta-pi);
%                 beta_bar1_dot=(sin(2*betaOA)/(0.5*((cos(2*betaOA))^2+1))-(Delta_beta_bar)/pi)*betaOA_dot;
                beta_bar1=beta_bar-pi;               
                end
                beta_bar1_dot=(b^(2*n)*(2*n-1)*cos(betaOA)*abs(cos(betaOA))^(2*n-2))/(a^(2*n)*sin(betaOA)*abs(sin(betaOA))^(2*n))*betaOA_dot;
            end
            F_AOk=  [cos(beta_bar1);sin(beta_bar1)];    
            F_A = F_A+ sigma*F_AOk;
            
            Px(i,j)=cos(beta_bar1);
            Py(i,j)=sin(beta_bar1);
       
        if E_AOk>EmAO
            dist=norm(rS-[X(i,j),Y(i,j)]');
            Px(i,j)=sigma*Px(i,j)+(1-sigma)*(rS(1)-X(i,j))/dist;
            Py(i,j)=sigma*Py(i,j)+(1-sigma)*(rS(2)-Y(i,j))/dist;
        end
    end
end
%%
%f1=@(x) [rS(1)-x(1)+,];
figure

xsup=a*(EmAO+1)^(1/(2*n))*cos(0:pi/50:pi/2).^(1/n);
ysup=b*(EmAO+1)^(1/(2*n))*sin(0:pi/50:pi/2).^(1/n);
xsup=[xsup,-flip(xsup),-xsup,flip(xsup)];
ysup=[ysup,flip(ysup),-ysup,-flip(ysup)];
xsup=xsup+rO(1,k);
ysup=ysup+rO(2,k);        
plot(xsup,ysup,'color',[0.5,0,1],'linewidth',0.5)
%set(ez,'color',[0.5,0,1])
hold on;
xsup=a*(EbarAO+1)^(1/(2*n))*cos(0:pi/50:pi/2).^(1/n);
ysup=b*(EbarAO+1)^(1/(2*n))*sin(0:pi/50:pi/2).^(1/n);
xsup=[xsup,-flip(xsup),-xsup,flip(xsup)];
ysup=[ysup,flip(ysup),-ysup,-flip(ysup)];
xsup=xsup+rO(1,k);
ysup=ysup+rO(2,k);        
plot(xsup,ysup,'color',[0.3,.91,0],'linewidth',0.5)
hold on;
xsup=a*(EuAO+1)^(1/(2*n))*cos(0:pi/50:pi/2).^(1/n);
ysup=b*(EuAO+1)^(1/(2*n))*sin(0:pi/50:pi/2).^(1/n);
xsup=[xsup,-flip(xsup),-xsup,flip(xsup)];
ysup=[ysup,flip(ysup),-ysup,-flip(ysup)];
xsup=xsup+rO(1,k);
ysup=ysup+rO(2,k);        
plot(xsup,ysup,'color',[0,.0,0],'linewidth',0.5)
hold on;
plot(rS(1),rS(2),'ro','markersize',12);
hold on;
legend(['$E_o=E_o^m=$',num2str(EmAO,3)],['$E_o=\bar E_o=$',num2str(EbarAO,3)],['$E_o=E_o^u=$',num2str(EuAO,3)],'$\textbf{r}_s$','location','northwest')
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

hf=quiver(X,Y,Px,Py,1.7)
set(hf,'color',[.8,.8,.8])
drawnow;
xlabel('x [m]')
ylabel('y [m]')
title('')
%text(rS(1)-1.2,rS(2)+0.95,['$r_s$(',num2str(rS(1)),',',num2str(rS(2)),')'],'fontsize',12,'backgroundcolor',[0.9,0.9,0.9])

if (0)
headSize=1.1;
fontSize=16;
arrowLength=3;
plot(rS(1),rS(2),'ro','markersize',8);
hold on;
text(rS(1)+0.52,rS(2)+0.15,['$\textbf{r}_s$(',num2str(rS(1)),', ',num2str(rS(2)),')'],'fontsize',12)
hold on;
%Tangent and field line at safe location
quiver(rS(1),rS(2),cos(betaOS)/1.2,sin(betaOS)/1.2,arrowLength,'linewidth',2,'color',[1 0 0],'MaxHeadSize',0.9*headSize)
quiver(rS(1),rS(2),cos(beta_barOS),sin(beta_barOS),arrowLength,'linewidth',2,'color',[1 0 .7],'MaxHeadSize',0.9*headSize)
rd=0.8;
dbeta_barOS=beta_barOS-betaOS;
plot(rS(1)+rd*cos(betaOS:dbeta_barOS/10:beta_barOS),rS(2)+rd*sin(betaOS:dbeta_barOS/10:beta_barOS),'color',[1 0 0.7])
rd=2.3;
text(rS(1)+rd*cos(beta_barOS/1.2),rS(2)+rd*sin(beta_barOS/1.2),'\boldmath$\Delta\bar\beta_{ok}$','fontsize',fontSize,'color',[1 0 0.7])
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


text(X(pt1,pt2)-2.02,Y(pt1,pt2)-.80,'$\textbf{r}(x,y)$','fontsize',12);
hold on;
%tangent
line([X(pt1,pt2),3+X(pt1,pt2)],[Y(pt1,pt2),Y(pt1,pt2)],'linestyle','--','color','b')
hold on;
quiver(X(pt1,pt2),Y(pt1,pt2),cos(Beta_bar(pt1,pt2)),sin(Beta_bar(pt1,pt2)),arrowLength,'linewidth',2,'color',[0.1 0.7 0.1],'MaxHeadSize',headSize)
hold on
rd=1.2;
plot(X(pt1,pt2)+rd*cos(0:Beta_bar(pt1,pt2)/10:Beta_bar(pt1,pt2)),Y(pt1,pt2)+rd*sin(0:Beta_bar(pt1,pt2)/10:Beta_bar(pt1,pt2)),'color',[0.1 0.7 0.1])
hold on;
rd=2.7;
text(X(pt1,pt2)+rd*cos(Beta_bar(pt1,pt2)/1.30),Y(pt1,pt2)+rd*sin(Beta_bar(pt1,pt2)/1.30),'\boldmath$\bar{{\beta}}_{ok}(\textbf{r})$','fontsize',fontSize,'color',[0.1 0.7 0.1])
%vector field
quiver(X(pt1,pt2),Y(pt1,pt2),Px(pt1,pt2),Py(pt1,pt2),arrowLength,'linewidth',2,'color',[1 0.5 0],'MaxHeadSize',headSize)
hold on;
rd=1.8;
plot(X(pt1,pt2)+rd*cos(0:Beta_bar1(pt1,pt2)/10:Beta_bar1(pt1,pt2)),Y(pt1,pt2)+rd*sin(0:Beta_bar1(pt1,pt2)/10:Beta_bar1(pt1,pt2)),'color',[1 0.5 0])
hold on
rd=2.2;
text(X(pt1,pt2)+rd*cos(Beta_bar1(pt1,pt2)/2.6),Y(pt1,pt2)+rd*sin(Beta_bar1(pt1,pt2)/2.6),'\boldmath$\bar\beta_{ok}^{\prime}(\textbf{r})$','fontsize',fontSize,'color',[1 0.5 0])
%current location
line([rO(1,1),7+rO(1,1)],[rO(2,1),rO(2,1)],'linestyle','--','color','b')
hold on
plot([rO(1,1),X(pt1,pt2)],[rO(2,1),Y(pt1,pt2)],'b')
hold on;
rd=1.15;
plot(rO(1,1)+rd*cos(0:BetaOF(pt1,pt2)/10:BetaOF(pt1,pt2)),rO(2,1)+rd*sin(0:BetaOF(pt1,pt2)/10:BetaOF(pt1,pt2)),'b')
hold on;
rd=1.5;
text(rO(1,1)+rd*cos(BetaOF(pt1,pt2)/5),rO(2,1)+rd*sin(BetaOF(pt1,pt2)/5),'\boldmath$\beta_{ok}(\textbf{r})$','fontsize',fontSize,'color','b')
%safe location
l=line([rO(1,1),rS(1)],[rO(2,1),rS(2)],'linestyle','--','color','r')
hold on
rd=4.5;
plot(rO(1,1)+rd*cos(0:betaOS/10:betaOS),rO(2,1)+rd*sin(0:betaOS/10:betaOS),'r')
hold on;
rd=4.8;
text(rO(1,1)+rd*cos(betaOS/2),rO(2,1)+rd*sin(betaOS/2),'\boldmath$\beta_{ok}(\textbf{r}_s)$','fontsize',fontSize,'color','r')
hold on
rd=2.65;
plot(rO(1,1)+rd*cos(BetaOF(pt1,pt2):(betaOS-BetaOF(pt1,pt2))/10:betaOS),rO(2,1)+rd*sin(BetaOF(pt1,pt2):(betaOS-BetaOF(pt1,pt2))/10:betaOS),'r')
hold on;
rd=3.3;
text(rO(1,1)+rd*cos(BetaOF(pt1,pt2)/1.1),rO(2,1)+rd*sin(BetaOF(pt1,pt2)/1.1),'\boldmath$\Delta\beta_{ok}$','fontsize',fontSize,'color','r')
hold on

xlabel('x [m]')
ylabel('y [m]')
title('')
 xlim([-1*xl,max(xl,1.1*rS(1))])
 ylim([-0.6*xl,1.3*rS(2)])
end
