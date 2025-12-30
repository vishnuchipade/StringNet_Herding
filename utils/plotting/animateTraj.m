function animateTraj(X0,rP,rS,rVO,rho_Acon,rho_S,t1,tlen,dt,phase_time,NA,ND,XD_Des0,WDString_mat)
% ------------------------------
% Author: Vishnu S. Chipade
% Email: vishnuc@umich.edu
% ------------------------------


%SetPlotDefaults;

circCos=cos(0:pi/50:2*pi);
circSin=sin(0:pi/50:2*pi);

global nO rho_P aO bO obs
rCO2=obs.rCO2;
N=NA+ND;%+1;
%tlen=length(time);
for ii=1:N
    X(2*ii-1:2*ii,:)=X0(4*ii-3:4*ii-2,:);  %get the positions only
end

XD=X(2*NA+1:end,:);


for j=1:ND
    XD_Des(2*j-1:2*j,:)=XD_Des0(4*j-3:4*j-2,:);
end
XAc=[sum(X(1:2:2*NA,:),1)/NA;sum(X(2:2:2*NA,:),1)/NA];

hfig=figure('units','normalized','outerposition',[.1 0.1 .7 .8]);
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset;
left = outerpos(1) + ti(1)+0.02;
bottom = outerpos(2) + ti(2)+0.04;
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4)-0.04;
ax.Position = [left bottom ax_width ax_height];

hold all;
markSize=1;
fontSize=18;
vid=VideoWriter('Results/test.avi');
vid.FrameRate=30;
open(vid);
makeVideo=1;


hold on;
plot(rP(1)+rho_P*cos(0:pi/100:2*pi), rP(2)+rho_P*sin(0:pi/100:2*pi),'--b','linewidth',0.5)
fill(rP(1)+rho_P*cos(0:pi/100:2*pi), rP(2)+rho_P*sin(0:pi/100:2*pi),'b','FaceAlpha',0.25)
hold on;
text(rP(1),rP(2),'$\mathcal{P}$','fontsize',fontSize);
hold on;
plot(rS(1)+rho_S*cos(0:pi/100:2*pi), rS(2)+rho_S*sin(0:pi/100:2*pi),'--g','linewidth',0.5)
fill(rS(1)+rho_S*cos(0:pi/100:2*pi), rS(2)+rho_S*sin(0:pi/100:2*pi),'g','FaceAlpha',0.25)
text(rS(1),rS(2),'$\mathcal{S}$','fontsize',fontSize);
%Intial positions and target and safe area
tic
for j=1:N
    if j<=NA
        plot(X(2*(j-1)+1,t1),X(2*(j-1)+2,t1),'ro','markersize',2*markSize);
        plot(X(2*(j-1)+1,t1),X(2*(j-1)+2,t1),'ro','markersize',8*markSize)
        hand(j)=plot(X(2*(j-1)+1,t1),X(2*(j-1)+2,t1),'ro','markersize',4*markSize);
    else
        plot(X(2*(j-1)+1,t1),X(2*(j-1)+2,t1),'bo','markersize',2*markSize);
        plot(X(2*(j-1)+1,t1),X(2*(j-1)+2,t1),'bo','markersize',8*markSize)
        if j<=NA+ND
            hand(j)=plot(X(2*(j-1)+1,t1),X(2*(j-1)+2,t1),'bo','markersize',4*markSize);
        else
            hand(j)=plot(X(2*(j-1)+1,t1),X(2*(j-1)+2,t1),'*','color',[0.5,0.5,.9],'markersize',4*markSize);
        end
    end
end
hold on;
handAc=plot(XAc(1,t1),XAc(2,t1),'mo','markersize',2*markSize);
handAcon=plot(XAc(1,t1)+rho_Acon*circCos,XAc(2,t1)+rho_Acon*circSin,'r--','markersize',2*markSize);
handTime=annotation('textbox', [0.9 .9 0 0], 'String', ['T = ',num2str(t1*dt),' s'], 'FitBoxToText', true);%text(1000,1000,['T = ',num2str(t1*dt),' s'])

if ND>0
    for j=1:ND
        handDDes(j)=plot(XD_Des(2*(j-1)+1,t1),XD_Des(2*(j-1)+2,t1),'bsquare','markersize',2*markSize);
        for jj=j+1:ND
            handDS(j,jj)=plot(XD([2*j-1,2*jj-1],t1)/WDString_mat(j,jj,t1),XD([2*j,2*jj],1)/WDString_mat(j,jj,t1),'b');
        end
    end
end
% if sum(WDString_mat(1,:,1))<2
%  handDS=plot(XD(1:2:end,1),XD(2:2:end,1));
% else
% handDS=plot([XD(1:2:end,1); XD(1,1)],[XD(2:2:end,1); XD(2,1)]);
% end
% flagStringNet=0;

global R_m_AO R_bar_AO R_u_AO R_m_AD R_bar_AD R_u_AD
global A_A_O B_A_O C_A_O D_A_O A_A_D B_A_D C_A_D D_A_D
global E_m_O E_bar_O E_u_O E_bar_AO E_u_AO
global E_m_DO E_bar_DO E_u_DO
%Plot obstacles
for k=1:size(rVO,2)
    %     w=rO(3,k);
    %     h=rO(4,k);
    
    %Effective area of barrier function
    E_u=E_u_DO;
    E_bar=E_bar_DO;
    E_m=E_m_DO;
    if(0)
        hold on;
        xsup=aO(k)*(E_u(k)+1)^(1/(2*nO(k)))*cos(0:pi/50:pi/2).^(1/nO(k));
        ysup=bO(k)*(E_u(k)+1)^(1/(2*nO(k)))*sin(0:pi/50:pi/2).^(1/nO(k));
        xsup=[xsup,-flip(xsup),-xsup,flip(xsup)];
        ysup=[ysup,flip(ysup),-ysup,-flip(ysup)];
        xsup=xsup+rCO2(1,k);
        ysup=ysup+rCO2(2,k);
        plot(xsup,ysup,'color',[0.3,1,.7])
        %Initial boundary for the combined vewctor field
        %         f2=@(x,y) 0.5*(2*(abs(x-rO(1,k)))/w)^(2*nO(k))+0.5*(2*(abs(y-rO(2,k)))/h)^(2*nO(k))-1-E_bar_O(k);
        %         ez2=ezplot(f2,[-xm+rO(1,k),xm+rO(1,k),-ym+rO(2,k),ym+rO(2,k)]);
        hold on;
        xsup=aO(k)*(E_bar(k)+1)^(1/(2*nO(k)))*cos(0:pi/50:pi/2).^(1/nO(k));
        ysup=bO(k)*(E_bar(k)+1)^(1/(2*nO(k)))*sin(0:pi/50:pi/2).^(1/nO(k));
        xsup=[xsup,-flip(xsup),-xsup,flip(xsup)];
        ysup=[ysup,flip(ysup),-ysup,-flip(ysup)];
        xsup=xsup+rCO2(1,k);
        ysup=ysup+rCO2(2,k);
        plot(xsup,ysup,'color',[0.3,1,.7])
        %Safe boundary for the formation
        %         f2=@(x,y) 0.5*(2*(abs(x-rO(1,k)))/w)^(2*nO(k))+0.5*(2*(abs(y-rO(2,k)))/h)^(2*nO(k))-1-E_m_O(k);
        %         ez2=ezplot(f2,[-xm+rO(1,k),xm+rO(1,k),-ym+rO(2,k),ym+rO(2,k)]);
        hold on;
        xsup=aO(k)*(E_m(k)+1)^(1/(2*nO(k)))*cos(0:pi/50:pi/2).^(1/nO(k));
        ysup=bO(k)*(E_m(k)+1)^(1/(2*nO(k)))*sin(0:pi/50:pi/2).^(1/nO(k));
        xsup=[xsup,-flip(xsup),-xsup,flip(xsup)];
        ysup=[ysup,flip(ysup),-ysup,-flip(ysup)];
        xsup=xsup+rCO2(1,k);
        ysup=ysup+rCO2(2,k);
        plot(xsup,ysup,'color',[0.3,1,.7])
    end
    if (0)
        hold on;
        xsup=aO(k)*(E_bar_AO(k)+1)^(1/(2*nO(k)))*cos(0:pi/50:pi/2).^(1/nO(k));
        ysup=bO(k)*(E_bar_AO(k)+1)^(1/(2*nO(k)))*sin(0:pi/50:pi/2).^(1/nO(k));
        xsup=[xsup,-flip(xsup),-xsup,flip(xsup)];
        ysup=[ysup,flip(ysup),-ysup,-flip(ysup)];
        xsup=xsup+rCO2(1,k);
        ysup=ysup+rCO2(2,k);
        plot(xsup,ysup,'color',[0.8,.1,.1])
        
        hold on;
        xsup=aO(k)*(E_u_AO(k)+1)^(1/(2*nO(k)))*cos(0:pi/50:pi/2).^(1/nO(k));
        ysup=bO(k)*(E_u_AO(k)+1)^(1/(2*nO(k)))*sin(0:pi/50:pi/2).^(1/nO(k));
        xsup=[xsup,-flip(xsup),-xsup,flip(xsup)];
        ysup=[ysup,flip(ysup),-ysup,-flip(ysup)];
        xsup=xsup+rCO2(1,k);
        ysup=ysup+rCO2(2,k);
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
    fill( XO{k},YO{k},[0.2 0.2 0.2])
    text(sum(rVO{k}(1,:))/length(rVO{k}(1,:)),sum(rVO{k}(2,:))/length(rVO{k}(2,:)),['$\mathcal{O}_',num2str(k),'$'],'fontsize',0.8*fontSize,'color',[1,1,1]);

   end
end
% for k=1:size(rO,2)
%     hold on;
%     plot(X2(k,:),Y2(k,:),'--k','linewidth',0.5)
% end


%plot m-air dimensions
%following variables denote the safe area considered inside the MAir
y_l = -36;
y_h = -2.5;
x_l = 0.2;
x_h = 22;
lx_mair = x_h-x_l;
ly_mair = y_h-y_l;
rect=rectangle('Position',[x_l,y_l,lx_mair,ly_mair],'LineStyle','--');
hold on
rect=rectangle('Position',[x_l-1,y_l-1,lx_mair+2,ly_mair+2]);
axis equal
xlabel('x');
ylabel('y');
xlim([-10,40]);
ylim([-50,5]);

title('');

if (1)
    
    
end


handText=text(0.7, .7, '','Units','normalized','fontsize',20);

handText2=text(0.4, .4, '','Units','normalized','fontsize',20);

FlagHistory=0;
ns=30;  %add or skip ns values everytime the plot is updated
dTi_phase=8/dt;
if FlagHistory==1
    for i=t1:ns:tlen
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
    for i=t1+1:ns:tlen
        for j=1:N
            if j<=NA
                set(hand(j),'XData',X(2*(j-1)+1,i),'YData',X(2*(j-1)+2,i));
            else
                set(hand(j),'XData',X(2*(j-1)+1,i),'YData',X(2*(j-1)+2,i));
            end
            hold on;
        end
        set(handAc,'XData',XAc(1,i),'YData',XAc(2,i))
        set(handAcon,'XData',XAc(1,i)+rho_Acon*circCos,'YData',XAc(2,i)+rho_Acon*circSin)
        countS=0;
        for j=1:ND
            %ind=find(WDString_mat(j,:,i)==1);
            set(handDDes(j),'XData',XD_Des(2*(j-1)+1,i),'YData',XD_Des(2*(j-1)+2,i));
            for jj=j+1:ND%length(ind)
                %set(handDS(j,ind(jj)),'XData',[XD(2*j-1,i),XD(2*ind(jj)-1,i)]/WDString_mat(j,ind(jj),i),'YData',[XD(2*j,i),XD(2*ind(jj),i)]/WDString_mat(j,ind(jj),i));
                set(handDS(j,jj),'XData',[XD(2*j-1,i),XD(2*jj-1,i)]/WDString_mat(j,jj,i),'YData',[XD(2*j,i),XD(2*jj,i)]/WDString_mat(j,jj,i));
                WDString_mat(j,jj,i)=0;
                WDString_mat(jj,j,i)=0;
            end
            
        end
        %Show time in window
        set(handTime,'String',['T = ',num2str(i*dt),' s'])
        
        if i>phase_time(1) && i<phase_time(1) + dTi_phase
            set(handText,'String', 'Gathering Phase is completed');%text(1000,1000,['T = ',num2str(t1*dt),' s'])
            
        elseif i>phase_time(2) && i<phase_time(2) + dTi_phase
            set(handText,'String', 'Seeking Phase is completed');
        elseif i>phase_time(3) && i<phase_time(3) + 2*dTi_phase
            set(handText,'String', 'Enclosing Phase is completed');
            set(handText2,'String','but 2 attackers escaped','Position',[0.7,.6,0])
        else
            set(handText,'String', '');
            set(handText2,'String', '');
        end
        
        if i>300/dt
            set(handText2,'String', 'The escaped attackers reached $\mathcal{P}$','Position',[0.33,.43,0]);
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
set(handText,'String', 'All enclosed attackers are herded to $\mathcal{S}$','Position',[0.64,.7,0]);

if makeVideo==1
    fr=getframe(gcf);
    writeVideo(vid,fr);
end
close(vid)
PlotTime=toc;
end
