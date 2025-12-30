function figNumber=plotTraj(X0,rP,rS,rVO,rho_S,time,NA,ND,markerStyle,lineStyle,figNumber,plotText,plotSafeArea,plotFinalPos)
% ------------------------------
% Author: Vishnu S. Chipade
% Email: vishnuc@umich.edu
% ------------------------------

%This functions considers rectangular obstacles
global REOmin REOmax REOcomb nO rho_P
global aO bO obs
rCO2=obs.rCO2;
plotAcm=0;
%%
xlim=50;
N=NA+ND;
tlen=length(time);
if figNumber==0
figNumber=figure('units','normalized','outerposition',[.15 0.2 .575 .75]);
axis equal
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2)+0.06;
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4)-0.06;
ax.Position = [left bottom ax_width ax_height];
else
    figure(figNumber)
end
markSize=1.5;
fontSize=18;
for ii=1:N
    X(2*ii-1:2*ii,:)=X0(4*ii-3:4*ii-2,:);  %get the positions only
end
XAc=[sum(X(1:2:2*NA,:),1)/NA;sum(X(2:2:2*NA,:),1)/NA];
tic
%Plot the protected area, safe area
if plotSafeArea
hold on;
plot(rP(1)+rho_P*cos(0:pi/100:2*pi), rP(2)+rho_P*sin(0:pi/100:2*pi),'--b','linewidth',0.5)
fill(rP(1)+rho_P*cos(0:pi/100:2*pi), rP(2)+rho_P*sin(0:pi/100:2*pi),'b','FaceAlpha',0.45)
hold on;
text(rP(1),rP(2),'$\mathcal{P}$','fontsize',fontSize);
hold on;
plot(rS(1)+rho_S*cos(0:pi/100:2*pi), rS(2)+rho_S*sin(0:pi/100:2*pi),'--g','linewidth',0.5)
fill(rS(1)+rho_S*cos(0:pi/100:2*pi), rS(2)+rho_S*sin(0:pi/100:2*pi),'g','FaceAlpha',0.45)
text(rS(1),rS(2),'$\mathcal{S}$','fontsize',fontSize);
end
%Intial positions and target
%dXD0=[-15,-40;-10,-45;-34,45;14,-25;14,15]';
%dXD0=[-15,-40;-10,-45;-34,-35;14,-25;14,15]';
dXD0=[15,-15;10,5;14,-15;14,5;14,15;-15,10]';
for j=1:N
    if j<=NA
        plot(X(2*(j-1)+1,1),X(2*(j-1)+2,1),'r.');
        hold on;
        hand(j)=plot(X(2*(j-1)+1,1),X(2*(j-1)+2,1),markerStyle,'color','r','markersize',5*markSize);
        %text(X(2*(j-1)+1,1),X(2*(j-1)+2,1),num2str(j))
    else
        plot(X(2*(j-1)+1,1),X(2*(j-1)+2,1),'b.');
        hold on;
        hand(j)=plot(X(2*(j-1)+1,1),X(2*(j-1)+2,1),markerStyle,'color','b','markersize',5*markSize);
        if plotText
        text(X(2*(j-1)+1,1)+dXD0(1,j-NA),X(2*(j-1)+2,1)+dXD0(2,j-NA),['$\mathcal{D}_',num2str(j-NA),'$'],'fontsize',fontSize)
        end
    end
    
end
if(plotAcm)
handAc=plot(XAc(1,1),XAc(2,1),'mo','markersize',8*markSize);
end
hold on;
global R_m_AO R_bar_AO R_u_AO R_m_AD R_bar_AD R_u_AD
global A_A_O B_A_O C_A_O D_A_O A_A_D B_A_D C_A_D D_A_D
global E_m_O E_bar_O E_u_O E_m_AO E_bar_AO E_u_AO E_v_DO
%Plot obstacles

for k=1:size(rVO,2)
    %     w=rO(3,k);
    %     h=rO(4,k);
    
    %Effective area of barrier function
    if(0)
        hold on;
        xsup=aO(k)*(E_u_O(k)+1)^(1/(2*nO(k)))*cos(0:pi/50:pi/2).^(1/nO(k));
        ysup=bO(k)*(E_u_O(k)+1)^(1/(2*nO(k)))*sin(0:pi/50:pi/2).^(1/nO(k));
        xsup=[xsup,-flip(xsup),-xsup,flip(xsup)];
        ysup=[ysup,flip(ysup),-ysup,-flip(ysup)];
        xsup=xsup+rO(1,k);
        ysup=ysup+rO(2,k);
        plot(xsup,ysup,'color',[0.3,1,.7])
        %Initial boundary for the combined vewctor field
        %         f2=@(x,y) 0.5*(2*(abs(x-rO(1,k)))/w)^(2*nO(k))+0.5*(2*(abs(y-rO(2,k)))/h)^(2*nO(k))-1-E_bar_O(k);
        %         ez2=ezplot(f2,[-xm+rO(1,k),xm+rO(1,k),-ym+rO(2,k),ym+rO(2,k)]);
        hold on;
        xsup=aO(k)*(E_bar_O(k)+1)^(1/(2*nO(k)))*cos(0:pi/50:pi/2).^(1/nO(k));
        ysup=bO(k)*(E_bar_O(k)+1)^(1/(2*nO(k)))*sin(0:pi/50:pi/2).^(1/nO(k));
        xsup=[xsup,-flip(xsup),-xsup,flip(xsup)];
        ysup=[ysup,flip(ysup),-ysup,-flip(ysup)];
        xsup=xsup+rO(1,k);
        ysup=ysup+rO(2,k);
        plot(xsup,ysup,'color',[0.3,1,.7])
        %Safe boundary for the formation
        %         f2=@(x,y) 0.5*(2*(abs(x-rO(1,k)))/w)^(2*nO(k))+0.5*(2*(abs(y-rO(2,k)))/h)^(2*nO(k))-1-E_m_O(k);
        %         ez2=ezplot(f2,[-xm+rO(1,k),xm+rO(1,k),-ym+rO(2,k),ym+rO(2,k)]);
        hold on;
        xsup=aO(k)*(E_m_O(k)+1)^(1/(2*nO(k)))*cos(0:pi/50:pi/2).^(1/nO(k));
        ysup=bO(k)*(E_m_O(k)+1)^(1/(2*nO(k)))*sin(0:pi/50:pi/2).^(1/nO(k));
        xsup=[xsup,-flip(xsup),-xsup,flip(xsup)];
        ysup=[ysup,flip(ysup),-ysup,-flip(ysup)];
        xsup=xsup+rO(1,k);
        ysup=ysup+rO(2,k);
        plot(xsup,ysup,'color',[0.3,1,.7])
        
        hold on;
        xsup=aO(k)*(E_v_DO(k)+1)^(1/(2*nO(k)))*cos(0:pi/50:pi/2).^(1/nO(k));
        ysup=bO(k)*(E_v_DO(k)+1)^(1/(2*nO(k)))*sin(0:pi/50:pi/2).^(1/nO(k));
        xsup=[xsup,-flip(xsup),-xsup,flip(xsup)];
        ysup=[ysup,flip(ysup),-ysup,-flip(ysup)];
        xsup=xsup+rO(1,k);
        ysup=ysup+rO(2,k);
        plot(xsup,ysup,'color',[0.1,1,.8])
        
    end
    
    if(0)
        hold on;
        xsup=aO(k)*(E_m_AO(k)+1)^(1/(2*nO(k)))*cos(0:pi/50:pi/2).^(1/nO(k));
        ysup=bO(k)*(E_m_AO(k)+1)^(1/(2*nO(k)))*sin(0:pi/50:pi/2).^(1/nO(k));
        xsup=[xsup,-flip(xsup),-xsup,flip(xsup)];
        ysup=[ysup,flip(ysup),-ysup,-flip(ysup)];
        xsup=xsup+rO(1,k);
        ysup=ysup+rO(2,k);
        plot(xsup,ysup,'color',[0.8,.1,.1])
        hold on;
    end
    
    if (0)
        hold on;
        xsup=aO(k)*(E_m_AO(k)+1)^(1/(2*nO(k)))*cos(0:pi/50:pi/2).^(1/nO(k));
        ysup=bO(k)*(E_m_AO(k)+1)^(1/(2*nO(k)))*sin(0:pi/50:pi/2).^(1/nO(k));
        xsup=[xsup,-flip(xsup),-xsup,flip(xsup)];
        ysup=[ysup,flip(ysup),-ysup,-flip(ysup)];
        xsup=xsup+rO(1,k);
        ysup=ysup+rO(2,k);
        plot(xsup,ysup,'color',[0.8,.1,.1])
        hold on;
        
        xsup=aO(k)*(E_bar_AO(k)+1)^(1/(2*nO(k)))*cos(0:pi/50:pi/2).^(1/nO(k));
        ysup=bO(k)*(E_bar_AO(k)+1)^(1/(2*nO(k)))*sin(0:pi/50:pi/2).^(1/nO(k));
        xsup=[xsup,-flip(xsup),-xsup,flip(xsup)];
        ysup=[ysup,flip(ysup),-ysup,-flip(ysup)];
        xsup=xsup+rO(1,k);
        ysup=ysup+rO(2,k);
        plot(xsup,ysup,'color',[0.8,.1,.1])
        
        hold on;
        xsup=aO(k)*(E_u_AO(k)+1)^(1/(2*nO(k)))*cos(0:pi/50:pi/2).^(1/nO(k));
        ysup=bO(k)*(E_u_AO(k)+1)^(1/(2*nO(k)))*sin(0:pi/50:pi/2).^(1/nO(k));
        xsup=[xsup,-flip(xsup),-xsup,flip(xsup)];
        ysup=[ysup,flip(ysup),-ysup,-flip(ysup)];
        xsup=xsup+rO(1,k);
        ysup=ysup+rO(2,k);
        plot(xsup,ysup,'color',[0.8,.1,.1])     
        
    end
    
    %For attacker
    if(0)
        f2=@(x,y) (abs(x-rO(1,k)))^(2)+(abs(y-rO(2,k)))^(2)-R_u_AO(k)^2;
        ez2=ezplot(f2,[-R_u_AO(k)+rO(1,k),R_u_AO(k)+rO(1,k),-R_u_AO(k)+rO(2,k),R_u_AO(k)+rO(2,k)]);
        set(ez2,'color','r')
        hold on;
    end
    %Obstacle (rectangular)
    %     w=rO(3,k);
    %     h=rO(4,k);
    if(0)
    XO{k}= [rVO{k}(1,:),rVO{k}(1,1)];
    YO{k}= [rVO{k}(2,:),rVO{k}(2,1)];
    hold on;
    fill(XO{k},YO{k},[0.2 0.2 0.2])
    text(rCO2(1,k)-7,rCO2(2,k),['$\mathcal{O}_',num2str(k),'$'],'fontsize',0.8*fontSize,'color',[1,1,1]);
    end
end
% for k=1:size(rO,2)
%     hold on;
%     plot(X2(k,:),Y2(k,:),'k','linewidth',0.5)
% end
xlabel('x [m]');
ylabel('y [m]');
title('');
%%
FlagHistory=1;
if FlagHistory==1
    for j=1:N
        if j<=NA
            plot(X(2*(j-1)+1,1:end),X(2*(j-1)+2,1:end),'r','markersize',markSize,'LineStyle',lineStyle)
%             plot(X(2*(j-1)+1,t_phase(1):t_phase(2)),X(2*(j-1)+2,t_phase(1):t_phase(2)),'r--','markersize',markSize)
%             plot(X(2*(j-1)+1,t_phase(2):t_phase(3)),X(2*(j-1)+2,t_phase(2):t_phase(3)),'ro','markersize',markSize)
%             plot(X(2*(j-1)+1,t_phase(3):end),X(2*(j-1)+2,t_phase(3):end),'r:','markersize',markSize)
            hold on;
            %quiver(X(2*(j-1)+1,1:end),X(2*(j-1)+2,1:end), FA_arr(2*(j-1)+1,:),FA_arr(2*(j-1)+2,:))
        else
            plot(X(2*(j-1)+1,1:end),X(2*(j-1)+2,1:end),'b--','markersize',markSize,'LineStyle',lineStyle)
%             plot(X(2*(j-1)+1,t_phase(1):t_phase(2)),X(2*(j-1)+2,t_phase(1):t_phase(2)),'b-','markersize',markSize)
%             plot(X(2*(j-1)+1,t_phase(2):t_phase(3)),X(2*(j-1)+2,t_phase(2):t_phase(3)),'b.-','markersize',markSize)
%             plot(X(2*(j-1)+1,t_phase(3):end),X(2*(j-1)+2,t_phase(3):end),'b.-','markersize',markSize)
        end
        hold on;
    end
     if(plotAcm)
    plot(XAc(1,:),XAc(2,:),'m--','markersize',markSize)
        end
end

drawnow;
%Final positions
if(plotFinalPos)
for j=1:N
    if j<=NA
        plot(X(2*(j-1)+1,tlen),X(2*(j-1)+2,tlen),'r^','markersize',5*markSize)
        text(X(2*(j-1)+1,tlen),X(2*(j-1)+2,tlen),num2str(j))
    else
        plot(X(2*(j-1)+1,tlen),X(2*(j-1)+2,tlen),'b^','markersize',5*markSize)
        text(X(2*(j-1)+1,tlen),X(2*(j-1)+2,tlen),num2str(j-NA))
    end
    hold on;
end
end
if(plotAcm)
plot(XAc(1,end),XAc(2,end),'md','markersize',markSize)
end
PlotTime=toc;

end
