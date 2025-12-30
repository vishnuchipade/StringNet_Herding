% ------------------------------
% Author: Vishnu S. Chipade
% Email: vishnuc@umich.edu
% ------------------------------

set(0,'defaultTextInterpreter','latex'); %trying to set the default
% ------------------------------
% Author: Vishnu S. Chipade
% Email: vishnuc@umich.edu
% ------------------------------


for ii=1:N
    X2(2*ii-1:2*ii,:)=X(4*ii-3:4*ii-2,:);  %get the positions only
end

XD=X2(2*NA+1:end,:);
markSize=8;
fontSize=18;

%Only the gathering phase
if (0)
    rVO=obs.rVO;
    if (1)
        fig=figure('units','normalized','outerposition',[.15 0.2 .575 .75]);
        fontSize=18;
        for k=1:size(rVO,2)
            rVO1{k}=[rVO{k},rVO{k}(:,1)];
            XO{k}= rVO1{k}(1,:);
            YO{k}= rVO1{k}(2,:);
            polyin=polyshape(rVO{k}(1,:),rVO{k}(2,:));
            [rCO2(1,k),rCO2(2,k)]=centroid(polyin);
            hold on;
            fill(XO{k},YO{k},[0.2 0.2 0.2])
            text(sum(rVO{k}(1,:))/length(rVO{k}(1,:)),sum(rVO{k}(2,:))/length(rVO{k}(2,:)),['$\mathcal{O}_',num2str(k),'$'],'fontsize',0.8*fontSize,'color',[1,1,1]);
        end
    end
    
    motionPlanForDefOpenForm(tanG,XD0,XD_des0,ND,1,fig.Number)
    savefig('Results/defender_goal_paths1.fig')
elseif(0)
    fig=openfig('Results/defender_goal_paths1.fig')
else
   % fig=figure('units','normalized','outerposition',[.15 0.2 .575 .75]);
end

%%
%Gathering + seeking phase

%gathering phase

%Plot the trajectories
fontSize=18;
ti1=1;
ti2=ti_g;
figNumber=plotTraj(X(:,ti1:ti2),rP,rS,rVO,rho_S,time(ti1:ti2),NA,ND,'o','--',0,1,1,0);
ti1=ti_g+1;
ti2=ti_s;
plotTraj(X(:,ti1:ti2),rP,rS,rVO,rho_S,time(ti1:ti2),NA,ND,'-','-',figNumber,0,0,0)
set(gca,'fontsize',fontSize)
hold on;

%Open formation in Gathering phase
for i=1:NA
    plot(X(4*i-3,ti_g),X(4*i-2,ti_g),'r*','markersize',markSize)
end
NA4=4*NA;
for j=1:ND
    plot(X(NA4+4*j-3,ti_g),X(NA4+4*j-2,ti_g),'b*','markersize',markSize)
    %line
    %plot([X(NA4+4*j-3,ti_1),X(NA4+4*j+1,ti_1)], [X(NA4+4*j+2,ti_1),X(NA4+4*j+2,ti_1)],'b*')
end
i=ti_g;
countS=0;
WDString_mat_temp=WDString_mat;
for j=1:ND
    ind=find(WDString_mat(j,:,i)==1);
    for jj=1:length(ind)
        countS=countS+1;
        plot([XD(2*j-1,i),XD(2*ind(jj)-1,i)]/WDString_mat_temp(j,ind(jj),i),[XD(2*j,i),XD(2*ind(jj),i)]/WDString_mat_temp(j,ind(jj),i),'b-','linewidth',2);
        plot([XD(2*j-1,i),XD(2*ind(jj)-1,i)]/WDString_mat_temp(j,ind(jj),i),[XD(2*j,i),XD(2*ind(jj),i)]/WDString_mat_temp(j,ind(jj),i),'--','color',[.6,.6,.6],'linewidth',1.2)
        WDString_mat_temp(j,ind(jj),i)=0;
        WDString_mat_temp(ind(jj),j,i)=0;
    end
end


%Open Formation after Seeking phase
if (1)
    for i=1:NA
        plot(X(4*i-3,ti_s),X(4*i-2,ti_s),'rd','markersize',markSize)
    end
    NA4=4*NA;
    for j=1:ND
        plot(X(NA4+4*j-3,ti_s),X(NA4+4*j-2,ti_s),'bd','markersize',markSize)
        %line
        %plot([X(NA4+4*j-3,ti_1),X(NA4+4*j+1,ti_1)], [X(NA4+4*j+2,ti_1),X(NA4+4*j+2,ti_1)],'b-')
    end
    i=ti_s+1;
    countS=0;
    WDString_mat_temp=WDString_mat;
    for j=1:ND
        ind=find(WDString_mat(j,:,i)==1);
        for jj=1:length(ind)
            countS=countS+1;
            plot([XD(2*j-1,i),XD(2*ind(jj)-1,i)]/WDString_mat_temp(j,ind(jj),i),[XD(2*j,i),XD(2*ind(jj),i)]/WDString_mat_temp(j,ind(jj),i),'b-','linewidth',2);
            plot([XD(2*j-1,i),XD(2*ind(jj)-1,i)]/WDString_mat_temp(j,ind(jj),i),[XD(2*j,i),XD(2*ind(jj),i)]/WDString_mat_temp(j,ind(jj),i),'--','color',[.6,.6,.6],'linewidth',1.2)
            WDString_mat_temp(j,ind(jj),i)=0;
            WDString_mat_temp(ind(jj),j,i)=0;
        end
    end
    
end

hold on
assignStr=[];
for i=1:length(assign)
    assignStr=[assignStr, num2str(assign(i)), ', '];
end
assignStr(end-1:end)=[];
% str={['$\circ \;$ Initial positions (@ t= ',num2str(0),' s)'],['$-- \;$ Paths during Gathering'],[ '$\ast \;$  Gathering (@ t= ',num2str(dt*ti_1),' s)'],...
%     ['$- \;$ Paths during Seeking'],['$\diamond \;$ Seeking (@ t= ',num2str(dt*ti_s),' s)']};
str={['$\circ \;$ Initial positions'],['$-- \;$ Paths during Gathering'],[ '$\ast \; $ Gathering (@ t= ',num2str(dt*ti_1),' s)'],['$- \;$ Paths during Seeking'],['$\diamond \;$ Seeking (@ t= ',num2str(dt*ti_s),' s)']};

annotation('textbox',[.13 .11 .29 .2],'string',str,'fontsize',fontSize,'interpreter','latex')

%%

%Enclosing +herding phase

% figNumber=plotTraj(X(:,ti1:ti2),rP,rS,rVO,rho_S,time(ti1:ti2),NA,ND,'-','-',0,1,1,0);
ti1=ti_s+1;
ti2=ti_e+45;
figNumber=plotTraj(X(:,ti1:ti2),rP,rS,rVO,rho_S,time(ti1:ti2),NA,ND,'-','-.',0,1,1,0)
ti1=ti_e+45;
ti2=length(time);
plotTraj(X(:,ti1:ti2),rP,rS,rVO,rho_S,time(ti1:ti2),NA,ND,'-',':',figNumber,0,0,1)
set(gca,'fontsize',fontSize)
hold on;


%Open Formation after Seeking phase
if (1)
    for i=1:NA
        plot(X(4*i-3,ti_s),X(4*i-2,ti_s),'rd','markersize',markSize)
    end
    NA4=4*NA;
    for j=1:ND
        plot(X(NA4+4*j-3,ti_s),X(NA4+4*j-2,ti_s),'bd','markersize',markSize)
        %line
        %plot([X(NA4+4*j-3,ti_1),X(NA4+4*j+1,ti_1)], [X(NA4+4*j+2,ti_1),X(NA4+4*j+2,ti_1)],'b-')
    end
    i=ti_s+1;
    countS=0;
    WDString_mat_temp=WDString_mat;
    for j=1:ND
        ind=find(WDString_mat(j,:,i)==1);
        for jj=1:length(ind)
            countS=countS+1;
            plot([XD(2*j-1,i),XD(2*ind(jj)-1,i)]/WDString_mat_temp(j,ind(jj),i),[XD(2*j,i),XD(2*ind(jj),i)]/WDString_mat_temp(j,ind(jj),i),'b-','linewidth',2);
            plot([XD(2*j-1,i),XD(2*ind(jj)-1,i)]/WDString_mat_temp(j,ind(jj),i),[XD(2*j,i),XD(2*ind(jj),i)]/WDString_mat_temp(j,ind(jj),i),'--','color',[.6,.6,.6],'linewidth',1.2)
            WDString_mat_temp(j,ind(jj),i)=0;
            WDString_mat_temp(ind(jj),j,i)=0;
        end
    end
    
end
%Closed Formation
if (1)
    for i=1:NA
        plot(X(4*i-3,ti1),X(4*i-2,ti1),'rp','markersize',markSize)
    end
    NA4=4*NA;
    for j=1:ND
        plot(X(NA4+4*j-3,ti1),X(NA4+4*j-2,ti1),'bp','markersize',markSize)
        %line
        %plot([X(NA4+4*j-3,ti_1),X(NA4+4*j+1,ti_1)], [X(NA4+4*j+2,ti_1),X(NA4+4*j+2,ti_1)],'b-')
    end
    i=ti_e+50;
    countS=0;
    WDString_mat_temp=WDString_mat;
    for j=1:ND
        ind=find(WDString_mat(j,:,i)==1);
        for jj=1:length(ind)
            countS=countS+1;
            plot([XD(2*j-1,i),XD(2*ind(jj)-1,i)]/WDString_mat_temp(j,ind(jj),i),[XD(2*j,i),XD(2*ind(jj),i)]/WDString_mat_temp(j,ind(jj),i),'b-','linewidth',2);
            plot([XD(2*j-1,i),XD(2*ind(jj)-1,i)]/WDString_mat_temp(j,ind(jj),i),[XD(2*j,i),XD(2*ind(jj),i)]/WDString_mat_temp(j,ind(jj),i),'--','color',[.6,.6,.6],'linewidth',1.2)
            WDString_mat_temp(j,ind(jj),i)=0;
            WDString_mat_temp(ind(jj),j,i)=0;
        end
    end
    
end
%herding
i=ti-2;
countS=0;
WDString_mat_temp=WDString_mat;
for j=1:ND
    ind=find(WDString_mat(j,:,i)==1);
    for jj=1:length(ind)
        countS=countS+1;
        plot([XD(2*j-1,i),XD(2*ind(jj)-1,i)]/WDString_mat_temp(j,ind(jj),i),[XD(2*j,i),XD(2*ind(jj),i)]/WDString_mat_temp(j,ind(jj),i),'b-','linewidth',2);
        plot([XD(2*j-1,i),XD(2*ind(jj)-1,i)]/WDString_mat_temp(j,ind(jj),i),[XD(2*j,i),XD(2*ind(jj),i)]/WDString_mat_temp(j,ind(jj),i),'--','color',[.6,.6,.6],'linewidth',1.2)
        WDString_mat_temp(j,ind(jj),i)=0;
        WDString_mat_temp(ind(jj),j,i)=0;
    end
end

hold on
assignStr=[];
for i=1:length(assign)
    assignStr=[assignStr, num2str(assign(i)), ', '];
end
assignStr(end-1:end)=[];
% str={['$\diamond \;$ Seeking (@ t= ',num2str(dt*ti_s),' s)'], ['$-. \;$ Paths during Enclosing'],['$\star \;$ StringNet (@ t= ',num2str(dt*ti_3),' s)']...
%     ['$: \;$ Paths during herding'],['$\triangle \;$  Final Positions (@ t= ',num2str(dt*ti),' s)']};
str={['$-. \;$ Paths during Enclosing'],['$\star \;$ StringNet (@ t= ',num2str(dt*ti_3),' s)'],['$.. \;$ Paths during herding'],['$\triangle \;$  Final Positions (@ t= ',num2str(dt*ti),' s)']};

annotation('textbox',[.56 .11 .345 .16],'string',str,'fontsize',fontSize,'interpreter','latex')

%%
% hold on
% assignStr=[];
% for i=1:length(assign)
%     assignStr=[assignStr, num2str(assign(i)), ', '];
% end
% assignStr(end-1:end)=[];
% str={['$\circ$ - Initial positions @ t= ',num2str(0),' s'],[ '$\ast$ - Gathering @ t= ',num2str(dt*ti_1),' s'],...
%     ['$\diamond$ - Seeking @ t= ',num2str(dt*ti_s),' s'], ['$\star$ - Closed StringNet @ t= ',num2str(dt*ti_3),' s'],...
%     ['$\triangle$ - Final Positions @ t= ',num2str(dt*ti),' s'], ['Assignment: [', assignStr, ']']};
%
% annotation('textbox',[.59 .12 .29 .205],'string',str,'fontsize',fontSize,'interpreter','latex')


