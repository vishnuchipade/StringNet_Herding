%This function plots the problem formulation for herding multiple attackers
% ------------------------------
% Author: Vishnu S. Chipade
% Email: vishnuc@umich.edu
% ------------------------------

%It plots obstacle populated environment, multiple attackers and defenders



set(0, 'defaultAxesTickLabelInterpreter','latex');
set(0, 'defaultLegendInterpreter','latex');
set(0, 'defaultTextInterpreter','latex');

%SetPlotDefaults;
fontSize=14;
%colors=distinguishable_colors(20);
colors={[0,0,1],[1,0,0],[0,1,0],[0,0,0],[0.4,1,0.1],[0.3,0.3,0.3],[0.5,0.9,0.5],[0,0.5,1],[.1,.51,1],[1,0.5,1],[0.5,0.5,1],[0.1,0.9,0.4],[0.9,0.4,0.5],[0,0.5,0.4]};

circCos=cos(0:pi/100:2*pi);
circSin=sin(0:pi/100:2*pi);

%Plot the figure for problem setup
rCO=[7,8;6,-13;-10,-10;-11,12]';
w=[4,7,3,6];
h=[4,6,6,4];
NO=length(w);
for k=1:NO
rVO{k}=rCO(:,k)+0.5*[-w(k),-h(k);w(k),-h(k);w(k),h(k);-w(k),h(k)]';
end
rVO{4}(:,1)=[];
rVO{2}(:,4)=[];
% rVO{4}(1,1)=0.5*(rVO{4}(1,1)+rVO{4}(1,2));
% rVO{4}(2,1)=0.5*(rVO{4}(2,1)+rVO{4}(2,4));
rho_P=4;
rS=[-1,18]';
rho_S=4;
rho_A=0.8;
rho_D=0.8;
rho_safe=1.9;
rho_A_s=4;
rA=[17,11;11,15;17,15;12,19]';
%rA=[16,10]';
vA=-rA;
NA=size(rA,2);
rD=[7,-6;0,-11;-17,9; -8,5;]';
vD=(rA-rD);
ND=size(rD,2);


markerSize=6;
headSize=2;
figure('units','normalized','outerposition',[.1 0.1 .7 .8])
axis equal
set(gca,'XTick',[],'YTick',[])
set(gca,'visible','off')
hold on 


%Plot defenders sensing zone
varRho_d=25;
fill(varRho_d*circCos,varRho_d*circSin,[ 0.5843 0.8157 0.9882])

varRho_da=22;
fill(varRho_da*circCos,varRho_da*circSin,[1,1,1])

%plot the protected area
plot(rho_P*circCos,rho_P*circSin,'color',colors{5})
fill(rho_P*circCos,rho_P*circSin,[0.8,0.8,1])
ang1=pi/4;
da=5*pi/180;
quiver(0,0,rho_P*cos(ang1),rho_P*sin(ang1),1.01,'color',colors{4})
text(0,-.4,['$\mathcal{P}$'],'fontsize',fontSize)
text(rho_P/2*cos(ang1-da),rho_P/2*sin(ang1-da),'$\rho_{pa}$','fontsize',fontSize)
text(-4.2,-1.8,'Safety-Critical','fontsize',0.8*fontSize)
text(-1.2,-2.8,'Area','fontsize',0.8*fontSize)

ang1=pi+pi/14;
da=5*pi/180;
plot([0,varRho_d*cos(ang1)],[0,varRho_d*sin(ang1)],'color',colors{4})
text(varRho_d/2*cos(ang1-da),varRho_d/2*sin(ang1-da),'$\varrho_d$','fontsize',fontSize)

ang1=pi-pi/20;
da=5*pi/180;
plot([0,varRho_da*cos(ang1)],[0,varRho_da*sin(ang1)],'color',colors{4})
text(varRho_da/2*cos(ang1-da),varRho_da/2*sin(ang1-da),'$\varrho_{d}^{a}$','fontsize',fontSize)



%plot the safe area
plot(rS(1)+rho_S*circCos,rS(2)+rho_S*circSin,'color',colors{5})
fill(rS(1)+rho_S*circCos,rS(2)+rho_S*circSin,colors{7})
ang1=pi/4;
da=5*pi/180;
quiver(rS(1),rS(2),rho_S*cos(ang1),rho_S*sin(ang1),1.01,'color',colors{4})
text(rS(1),rS(2)-.4,['$\mathcal{S}$'],'fontsize',fontSize)
text(rS(1)+rho_S/2*cos(ang1-da),rS(2)+rho_S/2*sin(ang1-da),'$\rho_{sa}$','fontsize',fontSize)
text(rS(1)-3.2,rS(2)-1.9,'Safe Area','fontsize',0.8*fontSize)

%Plot the obstacles
NO=length(rVO);
for k=1:NO
    nVO(k)=size(rVO{k},2);
    posVO=rVO{k};
    rCO(:,k)=1/nVO(k)*[sum(posVO,2)];
    posVO=[posVO,posVO(:,1)];
    rVO{k}=posVO;
    plot(posVO(1,:),posVO(2,:),'color',colors{3})
    
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
        AngVOVck{k}(ii)=atan2(yV-rCO(2,k),xV-rCO(1,k));
        if  AngVOVck{k}(ii)<0 
            AngVOVck{k}(ii)= AngVOVck{k}(ii)+2*pi;
        end
      
        if ii==1
              AngVOVck1(k)= AngVOVck{k}(ii);
        end
        AngVOVck{k}(ii)=AngVOVck{k}(ii)-AngVOVck1(k);
        if  AngVOVck{k}(ii)<0 
            AngVOVck{k}(ii)= AngVOVck{k}(ii)+2*pi;
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
        plot(xV+rho_safe*cos(ang1:(ang2-ang1)/50:ang2),yV+rho_safe*sin(ang1:(ang2-ang1)/50:ang2),'color',[0.7,0.7,0.7])
        plot(xV+rho_safe*cos(ang1:(ang2-ang1)/50:ang2),yV+rho_safe*sin(ang1:(ang2-ang1)/50:ang2),'--','color',colors{9})
        plot([xV,xV+rho_safe*cos(ang1)],[yV,yV+rho_safe*sin(ang1)],'--','color',[0.7,0.7,0.7])
        plot([xV,xV+rho_safe*cos(ang2)],[yV,yV+rho_safe*sin(ang2)],'--','color',[0.7,0.7,0.7])
        %plot the straight line
        plot(posVO2(1,2*ii:2*ii+1),posVO2(2,2*ii:2*ii+1),'color',colors{9})        
    end
    
end

for m=1:length(rVO)
    nV(m)=size(rVO{m},2);
    tempV=rVO{m};
    %tempV=[tempV,tempV(:,1)];
    rVO{m}=tempV;
    plot(tempV(1,:),tempV(2,:),'color',colors{4})
    fill(tempV(1,:),tempV(2,:),colors{6})
    %rcV=sum(tempV,1)/nV(m);
    text(rCO(1,m)-.5,rCO(2,m),['$\mathcal{O}_',num2str(m),'$'],'color',[1,1,1],'fontsize',fontSize)
end

%plot the attackers
drA=[0.4,-1;-3,-1;0.8,0;1,0.7]';
for i=1:NA
    plot(rA(1,i)+rho_A*circCos,rA(2,i)+rho_A*circSin,'color',colors{2});
    fill(rA(1,i)+rho_A*circCos,rA(2,i)+rho_A*circSin,colors{2});
    thetaA=atan2(vA(2,i),vA(1,i));
    rHead=repmat(rA(:,i),1,4)+0.8*rho_A*[cos(thetaA),sin(thetaA); -sin(thetaA),cos(thetaA); sin(thetaA),-cos(thetaA);cos(thetaA),sin(thetaA)]';
    plot(rHead(1,:),rHead(2,:),'color',colors{2})  
    fill(rHead(1,:),rHead(2,:),[1,1,1])  
    text(rA(1,i)+drA(1,i),rA(2,i)+drA(2,i),['$\mathcal{A}_',num2str(i),'$'],'fontsize',fontSize)
end
if NA~=1
text(11,23,'Attackers','color',[1,0,0],'fontsize',fontSize)
else
    text(12,16,'Attacker','color',[1,0,0],'fontsize',fontSize)
end
%plot the connectivity region around the attackers
rho_ac=7;
rAc=sum(rA,2)/NA;
plot(rAc(1)+rho_ac*circCos,rAc(2)+rho_ac*circSin,'r-.')
theta=pi/1.15;
plot([rAc(1),rAc(1)+rho_ac*cos(theta)],[rAc(2),rAc(2)+rho_ac*sin(theta)],'r-.')
text(rAc(1)+rho_ac/1.2*cos(theta),rAc(2)+rho_ac/1.2*sin(theta)+0.7,'$\rho_{ac}$','fontsize',fontSize)


%plot the defenders
for j=1:ND
    plot(rD(1,j)+rho_D*circCos,rD(2,j)+rho_D*circSin,'color',colors{1});
    fill(rD(1,j)+rho_D*circCos,rD(2,j)+rho_D*circSin,colors{1});
    thetaD=atan2(vD(2,j),vD(1,j));
    rHead=repmat(rD(:,j),1,4)+0.8*rho_D*[cos(thetaD),sin(thetaD); -sin(thetaD),cos(thetaD); sin(thetaD),-cos(thetaD);cos(thetaD),sin(thetaD)]';
    plot(rHead(1,:),rHead(2,:),'color',colors{1}) 
    fill(rHead(1,:),rHead(2,:),[1,1,1])  
    text(rD(1,j)+0.8,rD(2,j),['$\mathcal{D}_',num2str(j),'$'],'fontsize',fontSize)
end
text(-1,-8,'Defenders','color',[0,0,1],'fontsize',fontSize)

%Plot the projection of one of the defenders on the boundary of one of the
%obstacles (delta agent)
%the projection of D1 on actual obstacle
j=3;
k=4;
angj=atan2(rD(2,j)-rCO(2,k),rD(1,j)-rCO(1,k));
if angj<0
    angj=angj+2*pi;
end
angj=angj-AngVOVck1(k);
if angj<0
    angj=angj+2*pi;
end
indD=find(AngVOVck{k}(:)<angj);
ind=indD(end);
if (1)
    r1=rVO{k}(:,ind);
    r2=rVO{k}(:,ind+1);
    dx=r2(1)-r1(1);
    dy=r2(2)-r1(2);
    m=dy/dx;  %slope of line joining the defenders
    c=r1(2)-m*r1(1);
    rDProj(1,1)=(m*rD(2,j)+rD(1,j)-m*c)/(1+m^2);
    rDProj(2,1)=m*rDProj(1,1)+c;
    if m<1e16
        lambdaDP=(rDProj(1,1)-r1(1))/dx;
    else
        lambdaDP=(rDProj(2,1)-r1(1))/dy;
    end
    if lambdaDP<0
        rDProj(:,1)=r1;
    elseif lambdaDP>1
        rDProj(:,1)=r1;
    end
    rTP=r1-r2;
    rTP=rTP/norm(rTP);
    if rTP'*vD(:,j)<0
        rTP=-rTP;
    end   
end
%Projection on the outer boundary (approximated obstacle)
d=norm(rDProj-rD(:,j));
rDProj=(1-rho_safe/d)*rDProj+rho_safe/d*rD(:,j);

%Tangent at the projection point
vDProj=rTP'*vD(:,j)*rTP;

%Projection
plot(rDProj(1),rDProj(2),'o','color',colors{1},'markersize',2*markerSize);
text(rDProj(1)-2.5,rDProj(2)-.5,['$\mathbf{r}_{\delta', num2str(j) ,',',num2str(k),'}$'],'fontsize',fontSize,'color',[0,0,1])

quiver(rDProj(1),rDProj(2),vDProj(1),vDProj(2),0.1,'linewidth',2.5,'color',colors{1},'MaxHeadSize',.5*headSize)
text(rDProj(1)-0.3,rDProj(2)-1.8,['$\mathbf{v}_{\delta', num2str(j) ,',',num2str(k),'}$'],'fontsize',fontSize,'color',[0,0,1])


%Plot the projection of one of the attackers on the boundary of one of the
%obstacles (delta agent)
%the projection of D1 on actual obstacle
i=2;
k=1;
angj=atan2(rA(2,i)-rCO(2,k),rA(1,i)-rCO(1,k));
if angj<0
    angj=angj+2*pi;
end
angj=angj-AngVOVck1(k);
if angj<0
    angj=angj+2*pi;
end
indD=find(AngVOVck{k}(:)<angj);
ind=indD(end);
if (1)
    r1=rVO{k}(:,ind);
    r2=rVO{k}(:,ind+1);
    dx=r2(1)-r1(1);
    dy=r2(2)-r1(2);
    m=dy/dx;  %slope of line joining the defenders
    c=r1(2)-m*r1(1);
    if m<1e16
    rAProj(1,1)=(m*rA(2,i)+rA(2,i)-m*c)/(1+m^2);
    rAProj(2,1)=m*rAProj(1,1)+c;
    else
        rAProj(1,1)=rVO{k}(1,ind);
        rAProj(2,1)=rA(2,i);
    end
    if m<1e16
        lambdaDP=(rAProj(1,1)-r1(1))/dx;
    else
        lambdaDP=(rAProj(2,1)-r1(2))/dy;
    end
    if lambdaDP<0
        rAProj(:,1)=r1;
    elseif lambdaDP>1
        rAProj(:,1)=r1;
    end
    rTP=r1-r2;
    rTP=rTP/norm(rTP);
    if rTP'*vA(:,i)<0
        rTP=-rTP;
    end   
end
%Projection on the outer boundary (approximated obstacle)
d=norm(rAProj-rA(:,i));
rAProj=(1-rho_safe/d)*rAProj+rho_safe/d*rA(:,i);

%Tangent at the projection point
vAProj=rTP'*vA(:,i)*rTP;

%Projection
plot(rAProj(1),rAProj(2),'o','color',colors{2},'markersize',2*markerSize);
text(rAProj(1)+1.3,rAProj(2),['$\mathbf{r}_{\beta',num2str(i),',',num2str(k),'}$'],'fontsize',fontSize,'color',[1,0,0])

quiver(rAProj(1),rAProj(2),vAProj(1),vAProj(2),0.2,'linewidth',2.5,'color',colors{2},'MaxHeadSize',.5*headSize)
text(rAProj(1)+0.9,rAProj(2)-1.8,['$\mathbf{v}_{\beta',num2str(i),',',num2str(k),'}$'],'fontsize',fontSize, 'color',[1,0,0])


%Plot the sensing area of the attackers
plot(rA(1,3)+rho_A_s*circCos,rA(2,3)+rho_A_s*circSin,'r--')
ang=pi/1.5;
plot([rA(1,3),rA(1,3)+rho_A_s*cos(ang)],[rA(2,3),rA(2,3)+rho_A_s*sin(ang)],'r--')
text(rA(1,3)+rho_A_s*cos(ang)/1.72+0.25,rA(2,3)+rho_A_s*sin(ang)/1.72+0.2,'$\varrho_{a3}$','fontsize',fontSize)

%text for the vertices
k=2;
text(rVO{k}(1,1)+.35,rVO{k}(2,1)-0.95,['$\mathbf{r}_{o',num2str(k),'}^1$'],'fontsize',fontSize)
text(rVO{k}(1,2)+.15,rVO{k}(2,2)+0.85,['$\mathbf{r}_{o',num2str(k),'}^2$'],'fontsize',fontSize)
text(rVO{k}(1,3)+.15,rVO{k}(2,3)-0.95,['$\mathbf{r}_{o',num2str(k),'}^3$'],'fontsize',fontSize)

k=3;
for kk=1:2*nVO(k)
    plot(rVO2{k}(1,kk),rVO2{k}(2,kk),'o','color',colors{9})
end
text(rVO2{k}(1,1)-2.35,rVO2{k}(2,1),['$\mathbf{r}_{\bar{o}',num2str(k),'}^1$'],'fontsize',fontSize,'color',colors{9})
text(rVO2{k}(1,2),rVO2{k}(2,2)-0.95,['$\mathbf{r}_{\bar{o}',num2str(k),'}^2$'],'fontsize',fontSize,'color',colors{9})
text(rVO2{k}(1,3),rVO2{k}(2,3)-0.95,['$\mathbf{r}_{\bar{o}',num2str(k),'}^3$'],'fontsize',fontSize,'color',colors{9})
text(rVO2{k}(1,4)+.8,rVO2{k}(2,4),['$\mathbf{r}_{\bar{o}',num2str(k),'}^4$'],'fontsize',fontSize,'color',colors{9})
text(rVO2{k}(1,5)+0.8,rVO2{k}(2,5),['$\mathbf{r}_{\bar{o}',num2str(k),'}^5$'],'fontsize',fontSize,'color',colors{9})
text(rVO2{k}(1,6),rVO2{k}(2,6)+0.95,['$\mathbf{r}_{\bar{o}',num2str(k),'}^6$'],'fontsize',fontSize,'color',colors{9})
text(rVO2{k}(1,7),rVO2{k}(2,7)+0.95,['$\mathbf{r}_{\bar{o}',num2str(k),'}^7$'],'fontsize',fontSize,'color',colors{9})
text(rVO2{k}(1,8)-2.35,rVO2{k}(2,8),['$\mathbf{r}_{\bar{o}',num2str(k),'}^8$'],'fontsize',fontSize,'color',colors{9})

%radius of expansion
k=4;
plot([rVO{k}(1,1),rVO2{k}(1,1)],[rVO{k}(2,1),rVO2{k}(2,1)],'--','color',[0.2,0.2,0.2])
text((rVO{k}(1,1)+rVO2{k}(1,1))/2-1.9,(rVO{k}(2,1)+rVO2{k}(2,1))/2+.8,'$\rho_{oa}$','fontsize',fontSize)

%Axes
quiver(12,-2,1,0,3,'linewidth',1.5,'color',[0,0,0],'MaxHeadSize',.5*headSize)
text(14.9,-3,'$\hat{\mathbf{i}}$','fontsize',fontSize)
quiver(12,-2,0,1,3,'linewidth',1.5,'color',[0,0,0],'MaxHeadSize',.5*headSize)
text(11,0.5,'$\hat{\mathbf{j}}$','fontsize',fontSize)
%Workspace
text(14,-6,'$\mathcal{W}$','fontsize',1.5*fontSize)