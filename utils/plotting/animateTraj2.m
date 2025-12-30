function animateTraj2(X0,rP,rS,rO,rho_Fmax,rho_B,rho_S,time,NA,ND,XD_Des,WDString_mat)
% ------------------------------
% Author: Vishnu S. Chipade
% Email: vishnuc@umich.edu
% ------------------------------

%Plots on the figure 2 where tangent graph and paths are also plotted

global nO rho_P aO bO
N=NA+ND;
tlen=length(time);
for ii=1:N
    X(2*ii-1:2*ii,:)=X0(4*ii-3:4*ii-2,:);  %get the positions only
end

XD=X(2*NA+1:end,:);


for j=1:ND
    XD_Des(2*j-1:2*j,:)=X0(4*j-3:4*j-2,:);
end
XAc=[sum(X(1:2:2*NA,:),1)/NA;sum(X(2:2:2*NA,:),1)/NA];

figure(2)%,'units','normalized','outerposition',[.1 0.1 .7 .8])
markSize=1;
vid=VideoWriter('Results/test.avi');
vid.FrameRate=60;
open(vid);
makeVideo=1;

%Intial positions and target and safe area
tic
for j=1:N
    if j<=NA
        plot(X(2*(j-1)+1,1),X(2*(j-1)+2,1),'ro');
        hold on;
        hand(j)=plot(X(2*(j-1)+1,1),X(2*(j-1)+2,1),'ro','markersize',2*markSize);
    else
        plot(X(2*(j-1)+1,1),X(2*(j-1)+2,1),'bo');
        hold on;
        hand(j)=plot(X(2*(j-1)+1,1),X(2*(j-1)+2,1),'bo','markersize',2*markSize);
    end
    hold on;
    plot(rP(1),rP(2),'kp','markersize',6*markSize)
    hold on;
    plot(rP(1)+rho_P*cos(0:pi/100:2*pi), rP(2)+rho_P*sin(0:pi/100:2*pi),'--m','linewidth',0.5)
    hold on;
    plot(rS(1),rS(2),'go','markersize',markSize);
    hold on;
    plot(rS(1)+rho_S*cos(0:pi/100:2*pi), rS(2)+rho_S*sin(0:pi/100:2*pi),'--g','linewidth',0.5)
end
hold on;
handAc=plot(XAc(1,1),XAc(2,1),'mo','markersize',2*markSize);

for j=1:ND
    handD(j)=plot(XD_Des(2*(j-1)+1,1),XD_Des(2*(j-1)+2,1),'rsquare','markersize',2*markSize);
end
for j=1:ND-1
    handDS(j)=plot(XD(2*j-1:2:2*j+1,1)/WDString_mat(j,j+1,1),XD(2*j:2:2*j+2,1)/WDString_mat(j,j+1,1),'b');
end
handDS(ND)=plot([XD(2*ND-1,1),XD(1,1)]/WDString_mat(1,ND,1),[XD(2*ND,1),XD(2,1)]/WDString_mat(1,ND,1),'b');
% if sum(WDString_mat(1,:,1))<2
%  handDS=plot(XD(1:2:end,1),XD(2:2:end,1));
% else
% handDS=plot([XD(1:2:end,1); XD(1,1)],[XD(2:2:end,1); XD(2,1)]);
% end
% flagStringNet=0;

global R_m_AO R_bar_AO R_u_AO R_m_AD R_bar_AD R_u_AD
global A_A_O B_A_O C_A_O D_A_O A_A_D B_A_D C_A_D D_A_D
global E_m_O E_bar_O E_u_O E_bar_AO E_u_AO
%Plot obstacles
if(0)
for k=1:size(rO,2)
    w=rO(3,k);
    h=rO(4,k);
    
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
    end
    if (0)
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
    XO(k,:)= [rO(1,k)-w/2,rO(1,k)+w/2,rO(1,k)+w/2,rO(1,k)-w/2,rO(1,k)-w/2];
    YO(k,:)= [rO(2,k)-h/2,rO(2,k)-h/2,rO(2,k)+h/2,rO(2,k)+h/2,rO(2,k)-h/2];
    hold on;
    fill(XO(k,:),YO(k,:),[0.2 0.2 0.2])
end
end
% for k=1:size(rO,2)
%     hold on;
%     plot(X2(k,:),Y2(k,:),'--k','linewidth',0.5)
% end
axis equal
xlabel('x');
ylabel('y');
title('');

FlagHistory=0;
ns=30;  %add or skip ns values everytime the plot is updated
if FlagHistory==1
    for i=1:ns:tlen
        if i+ns-1<=tlen
            for j=1:N
                if j<=NA
                    plot(X(2*(j-1)+1,i:i+ns-1),X(2*(j-1)+2,i:i+ns-1),'ro','markersize',markSize);
                else
                    plot(X(2*(j-1)+1,i:i+ns-1),X(2*(j-1)+2,i:i+ns-1),'bo','markersize',markSize);
                end
                hold on;
            end
            plot(XAc(1,i:i+ns-1),XAc(2,i:i+ns-1),'mo','markersize',markSize)
            for j=1:ND
                plot(XD_Des(2*(j-1)+1,i:i+ns-1),XD_Des(2*(j-1)+2,i:i+ns-1),'rsquare','markersize',2*markSize);
            end
            
            drawnow;
            if makeVideo==1
                fr=getframe(gcf);
                writeVideo(vid,fr);
            end
        end
    end
elseif FlagHistory==0
    for i=2:ns:tlen
        for j=1:N
            if j<=NA
                set(hand(j),'XData',X(2*(j-1)+1,i),'YData',X(2*(j-1)+2,i));
            else
                set(hand(j),'XData',X(2*(j-1)+1,i),'YData',X(2*(j-1)+2,i));
            end
            hold on;
        end
        set(handAc,'XData',XAc(1,i),'YData',XAc(2,i))
        
        countS=0;
        for j=1:ND
            ind=find(WDString_mat(j,:,i)==1);
            
            for jj=1:length(ind)
                countS=countS+1;
            set(handDS(countS),'XData',[XD(2*j-1,i),XD(2*ind(jj)-1,i)]/WDString_mat(j,ind(jj),i),'YData',[XD(2*j,i),XD(2*ind(jj),i)]/WDString_mat(j,ind(jj),i));
            WDString_mat(j,ind(jj),i)=0;
            WDString_mat(ind(jj),j,i)=0;
            end
           
        end
       
        %         if sum(WDString_mat(1,:,i))<2 && flagStringNet~=1
        %             set(handDS,'XData',XD(1:2:end,i),'YData', XD(2:2:end,i));
        %         else
        %             if flagStringNet~=1
        %                 delete(handDS);
        %                 clear handDS;
        %                 handDS=plot([XD(1:2:end,i); XD(1,i)],[XD(2:2:end,i); XD(2,i)]);
        %                 flagStringNet=1;
        %             end
        %             set(handDS,'XData',[XD(1:2:end,i); XD(1,i)],'YData', [XD(2:2:end,i); XD(2,i)]);
        %         end
        drawnow;
        if makeVideo==1
            fr=getframe(gcf);
            writeVideo(vid,fr);
        end
        %drawnow;
    end
end

%Final positions
for j=1:N
    if j<=NA
        plot(X(2*(j-1)+1,i:tlen),X(2*(j-1)+2,i:tlen),'ro','markersize',markSize);
        hold on;
        plot(X(2*(j-1)+1,tlen),X(2*(j-1)+2,tlen),'rd','markersize',8)
    else
        plot(X(2*(j-1)+1,i:tlen),X(2*(j-1)+2,i:tlen),'bo','markersize',markSize);
        hold on;
        plot(X(2*(j-1)+1,tlen),X(2*(j-1)+2,tlen),'bd','markersize',8)
    end
    hold on;
    plot(XAc(1,i:tlen),XAc(2,i:tlen),'mo','markersize',markSize);
    hold on;
    plot(XAc(1,tlen),XAc(2,tlen),'md','markersize',markSize)
end
if makeVideo==1
    fr=getframe(gcf);
    writeVideo(vid,fr);
end
close(vid)
PlotTime=toc;
end
