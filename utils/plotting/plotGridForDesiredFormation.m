% ------------------------------
% Author: Vishnu S. Chipade
% Email: vishnuc@umich.edu
% ------------------------------

%SetPlotDefaults;
set(groot,'defaulttextinterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultTextInterpreter','latex');
set(groot, 'DefaultLineLineWidth', 2);
set(groot, 'DefaultAxesFontSize', 18);

colors={[0,0,1],[1,0,0],[0,0,0],[0.5,0.5,0.5],[0,0.5,1],[1,0.5,1],[0.5,0.5,1],[1,.1,0.5],[0.1,0.9,0.4],[0.9,0.4,0.5],[0,0.5,0.4],[1,0.4,0.4]};
fontSize=20;

rVO={[1,1;8,1;8,6;1,6]',[-3,17;3,17;3,23;-3,23;]'};
rV0{1}=[rVO{1}, rVO{1}(:,1)];
rV0{2}=[rVO{2}, rVO{2}(:,1)];

rCO(1,1)=(rVO{1}(1,1)+rVO{1}(1,2))/2;
rCO(2,1)=(rVO{1}(2,2)+rVO{1}(2,3))/2;
global rho_safe rho_D
rho_safe=1;
rho_D=0.2;
NO=length(rVO);
%find the tangent graph
tanG=tangentGraph(rVO);
%%

fig=figure('units','normalized','outerposition',[.2 0.2 .5 .8]);
hold all;
x=[-8.5:1.5:17];
y=[-8.5:1.5:17];
[xgrid,ygrid]=meshgrid(x,y);
xgrid=xgrid(:);
ygrid=ygrid(:);
ind1=find(xgrid<0);
ind2=find(xgrid>8);

ind3=find(ygrid<0);
ind4=find(ygrid>7);
ind=[ind1;ind2;ind3;ind4];

plot(xgrid(ind),ygrid(ind),'k.')

xlim([-6.2,17]);
ylim([-8.4,17]);

plot(rV0{1}(1,:),rV0{1}(2,:),'k-')
hold all
for k=1:NO-1
    nVO(k)=size(rVO{k},2);
    posVO=rVO{k};
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
            
            countV=countV+1;
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
            
            countV=countV+1;
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
        %corresponding to the next vertex
        posVO2(:,countV)=[xip,yip]';
    end
    %Shift the outer vertex corresponding to the first inner vertex to
    %appropriate position
    posVO2=posVO2(:,[countV,1:countV-1]);
    rVO2{k}=posVO2;
    nV2=countV;
    posVO2=[posVO2,posVO2(:,1)];
    
    %Find the perimeter and Plot the outer vertices
    PeriO(k)=0;
    for ii=1:nVO(k)
        %Get the inner vertex coordinates
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
        
        %Find the perimeters of the segments of the approximated obstacle
        PeriSOk{k}(2*ii-1,1)=rho_safe*(ang2-ang1);
        PeriSOk{k}(2*ii,1)=norm(posVO2(:,2*ii)-posVO2(:,2*ii+1));
        
        AngOk{k}(2*ii-1,1)=ang1;
        AngOk{k}(2*ii,1)=ang2;
        
        %plot the circular arc
        plot(xV+rho_safe*cos(ang1:(ang2-ang1)/50:ang2),yV+rho_safe*sin(ang1:(ang2-ang1)/50:ang2),'color',colors{4})
        %plot the straight line
        plot(posVO2(1,2*ii:2*ii+1),posVO2(2,2*ii:2*ii+1),'color',colors{4})
        
    end
end


%%

ND=3;
rI=[-0.1,-3.3;-4.1,-.1;8.75,-7 ]';
rF0=[13,13]';
RF=3.7;
theta0=pi/2.9;
for i=1:ND
    theta=theta0+pi/2+(i-1)*pi/(ND-1);
    rF(:,i)=rF0+RF*[cos(theta),sin(theta)]';
end

rIgrid=[0.5,-4;-4,0.5;9.5,-7]';
rFgrid=[9.5,14;11,9.5;15.5,11]';

%plot the points
markerSize=6;

colors={[0.5,0.5,1],[0,0.5,1],[0.5,0,1],[1,0.1,0.4],[.7,0.2,0.3]};
%Find the shortest paths and plot
for i=1:ND
    for ii=1:ND
        [Path(i,ii),tanG_prime1(i,ii),rE_dash1] = findPath(tanG,rIgrid(:,i),rFgrid(:,ii),0,[]);
        PlotShortestPath(tanG,Path(i,ii),tanG_prime1,fig.Number,[colors{i},1],'-',1)
    end
end
i1=1;
ii1=3;
i2=3;
ii2=2;
plotPathIntersections(Path(i1,ii1),Path(i2,ii2),tanG,tanG_prime1(i1,ii1),tanG_prime1(i2,ii2),fig.Number,{[colors{i1},0.2],[colors{i2},0.2]},12)
dx=rFgrid(1,ii2)-rIgrid(1,i2);
dy=rFgrid(2,ii2)-rIgrid(2,i2);
dx=dx/norm([dx,dy]);
dy=dy/norm([dx,dy]);
kappa=0.06;
annotation('doublearrow','Position',[0.71,0.44,kappa*dx,kappa*dy],'color',colors{i2},'linewidth',1.5);
text(11.5,3,'$\mathcal{L}_{3,2}^{1,3}$','color',colors{i2},'fontsize',fontSize)

%find the circle tangent to the shortest path and passing through the
%initial point
%%
if(0)
init=[2,1;1,1;1,1;]';
for i=1:ND
    m1=(Path(i,1).rV(2,2)-Path(i,1).rV(2,1))/(Path(i,1).rV(1,2)-Path(i,1).rV(1,1));
    a1=m1;
    b1=-1;
    d1=-m1*Path(i,1).rV(1,1) +Path(i,1).rV(2,1);
    %r1=abs(a1*rI(1,i)+b1*rI(2,i)+d1)/sqrt(a1^2+b1^2);
    f=@(x) [(rI(1,i)-x(1))^2+(rI(2,i)-x(2))^2-x(3)^2; (Path(i,1).rV(1,1)-x(1))^2+(Path(i,1).rV(2,1)-x(2))^2-x(3)^2;  x(3)-abs(a1*x(1)+b1*x(2)+d1)/sqrt(a1^2+b1^2)];
    rc1=fsolve(f,[rI(:,i);0.5]+[init(:,i);0]);
    ang1=atan2(rI(2,i)-rc1(2),rI(1,i)-rc1(1));
    if ang1<0
        ang1=ang1+2*pi;
    end
    ang2=atan2(Path(i,1).rV(2,1)-rc1(2),Path(i,1).rV(1,1)-rc1(1));
    if ang2<0
        ang2=ang2+2*pi;
    end
    ang=ang1:(ang2-ang1)/50:ang2;
    
    plot(rc1(1,1)+rc1(3)*cos(ang),rc1(2,1)+rc1(3)*sin(ang),'--','color',colors{i})
end
end


for i=1:ND
    plot([rI(1,i),rIgrid(1,i)],[rI(2,i),rIgrid(2,i)],'--','color',colors{i})
    plot(rI(1,i),rI(2,i),'bo','markersize',markerSize);
    plot(rIgrid(1,i),rIgrid(2,i),'d','markersize',markerSize,'color',colors{4});
    plot(rF(1,i),rF(2,i),'kp','markersize',markerSize);
    plot(rFgrid(1,i),rFgrid(2,i),'square','markersize',markerSize,'color',colors{5});
end
theta_arr=theta0+pi/2:pi/100:theta0+3*pi/2;
plot(rF0(1)+RF*cos(theta_arr), rF0(1)+RF*sin(theta_arr),'k--')

dr=RF/5*[cos([theta0+pi/2,theta0+3*pi/2]);sin([theta0+pi/2,theta0+3*pi/2])]
plot(rF0(1)+dr(1,:),rF0(2)+dr(2,:),'k-')
plot([rF0(1),rF0(1)+2.5],[rF0(2),rF0(2)],'k--')
theta_arr=0:theta0/50:theta0;
dR=1.3;
plot(rF0(1)+dR*cos(theta_arr), rF0(1)+dR*sin(theta_arr),'k--')
text(rF0(1)+dR,rF0(2)+.6*dR,'$\phi$','fontsize',fontSize)
quiver(rF0(1),rF0(2),cos(theta0),sin(theta0),3,'color',[0,0,0],'maxheadsize',2)
text(rF0(1)-.5,rF0(2)-.6,'$\mathbf{r}_{df}$','fontsize',fontSize)

%%
hold on
drI=[-0.5,1.3;-0.5,-1.3;-0.5,-1.3;]';
for i=1:ND
text(rI(1,i)+drI(1,i),rI(2,i)+drI(2,i),['$\mathbf{r}_{d',num2str(i),'}$'],'fontsize',fontSize,'color',[0,0,1])
end

drIgrid=[-0.5,-1.3;-0.5,1.3;0.5,-1;]';
for i=1:ND
text(rIgrid(1,i)+drIgrid(1,i),rIgrid(2,i)+drIgrid(2,i),['$\check{\mathbf{r}}_{d',num2str(i),'}$'],'fontsize',fontSize,'color',[1,0.1,0.4])
end

drF=[-1,.3;0.27,1;0.3,-.8;]';
for i=1:ND
text(rF(1,i)+drF(1,i),rF(2,i)+drF(2,i),['$\mathbf{\xi}_{',num2str(i),'}^{sc}$'],'fontsize',fontSize,'color',[0,0,0])
end

drFgrid=[.3,-.2;-1.0,-.55;-0.7,.9;]';
for i=1:ND
text(rFgrid(1,i)+drFgrid(1,i),rFgrid(2,i)+drFgrid(2,i),['$\mathbf{\xi}_{',num2str(i),'}^g$'],'fontsize',fontSize,'color',colors{5})
end



%Line joining initial position and the center of the obstacle
m1=(rCO(2,1)-rI(2,1))/(rCO(1,1)-rI(1,1));
c1=rI(2,1)-m1*rI(1,1);

x1=-6:0.2:15;
y1=m1.*x1+c1;

% plot(x1,y1,'-.','color',[0.5,0.5,0.9])
% plot(rCO(1),rCO(2),'*')



function [Path,tanG_prime,rE_dash] = findPath(tanG,rI,rF,flag,path0)

%This function evaluates a shortest path between the points rI and rF and
%static graph G0

%Input: tangG={G0, G0_pathType, rVO, rVO2, gTO_all, rTO_all, obsId_all,
%tangentIdO_all}
%Output:
global rho_safe
options=optimoptions('fmincon','display','off','algorithm','interior-point');%sqp-legacy');
colors={[0,0,1],[1,0,0],[0,0,0],[0.5,0.5,0.5],[0,0.5,1],[1,0.5,1],[0.5,0.5,1],[1,.1,0.5],[0.1,0.9,0.4],[0.3,0.4,0.5],[1,0.1,0.1]};

gammaMax=0.99999999999999;

NO=tanG.NO;
nVO=tanG.nVO;
xOk=tanG.xOk;
yOk=tanG.yOk;
dxOk_dg=tanG.dxOk_dg;
dyOk_dg=tanG.dyOk_dg;
GammaOk=tanG.GammaOk;
PeriSOk=tanG.PeriSOk;
PeriO=tanG.PeriO;
G=tanG.G;
G_pathType=tanG.G_pathType;
rVO=tanG.rVO;
rVO2=tanG.rVO2;
gTO_all=tanG.gTO_all;
rTO_all=tanG.rTO_all;
obsId_all=tanG.obsId_all;
obsVertId_all= tanG.obsVertId_all;
obsVertPos_all= tanG.obsVertPos_all;



%Ellipse filter
rE_dash=[];
if(0)
    E_m=tanG.E_m;
    %standard ellipse parameters
    cE=norm(rF-rI)/2;
    aE=E_m*cE;
    bE=cE*sqrt(E_m^2-1);
    xE=aE*cos(0:pi/100:2*pi);
    yE=bE*sin(0:pi/100:2*pi);
    rCE=(rI+rF)/2;
    angE=atan2(rF(2)-rI(2),rF(1)-rI(1));
    rotMatE=[cos(angE), -sin(angE);  sin(angE), cos(angE)];
    [rE_dash]=rotMatE*[xE;yE];
    rE_dash=rE_dash+rCE;
    
end

NTV0=size(G,2);  %Initial number of tangent vertices
indIF(1)=NTV0+1;
indIF(2)=NTV0+2;

%Check if the line between the initial and final point interesect any
%obstacles. If not then add an edge between them
r1=rI;
r2=rF;
dr=r2-r1;
drx=r2(1)-r1(1);
dry=r2(2)-r1(2);
mL1=(r2(2)-r1(2))/(r2(1)-r1(1));
cL1=r1(2)-mL1*r1(1);
dr_hat=dr/norm(dr);
flagRemoveIF=0; %for line between initial and final position
for k=1:NO
    %if k~=ki1 && k~=ki2
    %countCr=0;
    %check with the circles at the inner vertices
    for i=1:nVO(k)
        P1V=rVO{k}(1:2,i)-r1;
        distLV=abs(dr_hat(1)*P1V(2)-dr_hat(2)*P1V(1));
        if distLV<rho_safe  %if lines intersects the circle at the vertex
            %tangent_rPO(kt)=[];
            b=-2*(rVO{k}(1,i)-mL1*cL1+mL1*rVO{k}(2,i));
            a=(1+mL1^2);
            c=rVO{k}(1,i)^2+cL1^2+rVO{k}(2,i)^2-2*cL1*rVO{k}(2,i)-rho_safe^2;
            x_int=(-b+sqrt(b^2-4*a*c))/(2*a);
            lambda1=(x_int-r1(1))/drx;
            if (0<= lambda1) && (lambda1<=1)  % %if line segment intersects the line joining two vertices
                flagRemoveIF=1;
                break;
            end
        end
    end
    if flagRemoveIF==1
        break;
    end
    %check with the lines between the inner vertices
    for i=1:nVO(k)
        mL2=(rVO{k}(2,i)-rVO{k}(2,i+1))/(rVO{k}(1,i)-rVO{k}(1,i+1));
        cL2=rVO{k}(2,i)-mL2*rVO{k}(1,i);
        drx2=rVO{k}(1,i+1)-rVO{k}(1,i);
        x_int=(cL2-cL1)/(mL1-mL2);
        %y_int=mL1*x_int+cL1;
        lambda1=(x_int-r1(1))/drx;
        lambda2=(x_int-rVO{k}(1,i))/drx2;
        if (0<= lambda1) && (lambda1<=1) && (0<= lambda2) && (lambda2<=1) % %if line segments intersects the line joining two vertices
            flagRemoveIF=1;
            break;
        end
    end
    %end
    if flagRemoveIF==1
        break;
    end
end

if flagRemoveIF==0
    G(indIF(1),indIF(2))=norm(rI-rF);
    G_pathType(indIF(1),indIF(2))=3;
    
end

%Now update the graph G0 by adding tangents from rI and rF to each of the
%obstacles and then apply the dijkstra on the updated graph



NTVP=0;
%find tangents from rI and rF to all the obstacles
A=[eye(NO);-eye(NO)];
for rp=1:2
    if rp==1
        rP=rI;
    elseif rp==2
        rP=rF;
    end
    countP=1;
    countTotT=0;
    for k=1:NO
        m=@(g) (yOk{k}(g) - rP(2))/(xOk{k}(g) - rP(1));
        m1=@(g) dyOk_dg{k}(g)/dxOk_dg{k}(g);
        f=@(g) (m(g)-m1(g))^2;
        
        countT=0;
        
        for i=1:nVO(k)
            b=zeros(2*NO,1);
            b(k,1)=[GammaOk{k}(2*i)]';
            b(NO+k,1)=[-GammaOk{k}(2*i-1)]';
            g0=zeros(1,NO);
            alp1=.9;
            g0(k)=alp1*GammaOk{k}(2*i-1)+(1-alp1)*GammaOk{k}(2*i);
            %g0(kk)=alp1*GammaOk{kk}(2*j-1)+(1-alp1)*GammaOk{kk}(2*j);
            gt0=fmincon(f,g0,A,b,[],[],[],[],[],options);
            ftol=1e-7;
            gtol=1e-5;
            % gt0=fmincon(f,[GammaOk{k}(2*i-1),GammaOk{kk}(2*j-1)],A,b);
            if f(gt0)<ftol && GammaOk{k}(2*i-1)<gt0(k) && gt0(k)<GammaOk{k}(2*i)
                if countT==0
                    countT=countT+1;
                    gt(:,countT)=gt0;
                elseif ~sum( (abs(gt(k,:)-gt0(k))<gtol))
                    countT=countT+1;
                    gt(:,countT)=gt0;
                end
            end
            alp2=0.20;
            g0(k)=alp2*GammaOk{k}(2*i-1)+(1-alp2)*GammaOk{k}(2*i);
            %g0(kk)=alp2*GammaOk{kk}(2*j-1)+(1-alp2)*GammaOk{kk}(2*j);
            gt0=fmincon(f,g0,A,b,[],[],[],[],[],options);
            % gt0=fmincon(f,[GammaOk{k}(2*i-1),GammaOk{kk}(2*j-1)],A,b);
            if f(gt0)<ftol && GammaOk{k}(2*i-1)<gt0(k) && gt0(k)<GammaOk{k}(2*i)
                if countT==0
                    countT=countT+1;
                    gt(:,countT)=gt0;
                elseif ~sum((abs(gt(k,:)-gt0(k))<gtol))
                    countT=countT+1;
                    gt(:,countT)=gt0;
                end
            end
            
            g0(k)=alp1*GammaOk{k}(2*i-1)+(1-alp1)*GammaOk{k}(2*i);
            % g0(kk)=alp2*GammaOk{kk}(2*j-1)+(1-alp2)*GammaOk{kk}(2*j);
            gt0=fmincon(f,g0,A,b,[],[],[],[],[],options);
            % gt0=fmincon(f,[GammaOk{k}(2*i-1),GammaOk{kk}(2*j-1)],A,b);
            if f(gt0)<ftol && GammaOk{k}(2*i-1)<gt0(k) && gt0(k)<GammaOk{k}(2*i)
                if countT==0
                    countT=countT+1;
                    gt(:,countT)=gt0;
                elseif ~sum( (abs(gt(k,:)-gt0(k))<gtol))
                    countT=countT+1;
                    gt(:,countT)=gt0;
                end
            end
            
            g0(k)=alp2*GammaOk{k}(2*i-1)+(1-alp2)*GammaOk{k}(2*i);
            %g0(kk)=alp1*GammaOk{kk}(2*j-1)+(1-alp1)*GammaOk{kk}(2*j);
            gt0=fmincon(f,g0,A,b,[],[],[],[],[],options);
            % gt0=fmincon(f,[GammaOk{k}(2*i-1),GammaOk{kk}(2*j-1)],A,b);
            if f(gt0)<ftol && GammaOk{k}(2*i-1)<gt0(k) && gt0(k)<GammaOk{k}(2*i)
                if countT==0
                    countT=countT+1;
                    gt(:,countT)=gt0;
                elseif ~sum( (abs(gt(k,:)-gt0(k))<gtol))
                    countT=countT+1;
                    gt(:,countT)=gt0;
                end
            end
            
        end
        
        %X=-2:0.1:15;
        for i=1:countT
            rTPO(:,i)=[xOk{k}([gt(:,i)']),yOk{k}(gt(:,i)')];
            %rF(:,i)=[xOk{kk}([gt(:,i)']),yOk{kk}(gt(:,i)')];
            
            % Store the tangent points on the curves
            %             countTO(k)= countTO(k)+1;
            %             countTO(kk)= countTO(kk)+1;
            %             rTO{k}(:,countTO(k))=rP(:,i);
            %             rTO{kk}(:,countTO(kk))=rF(:,i);
            %             gTO{k}(countTO(k))=gt(k,i);
            %             gTO{kk}(countTO(kk))=gt(kk,i);
            
            %Store the information of tangents globally
            tangent_rPO(1:2,countTotT+i)=rTPO(:,i);
            tangent_rPO(3,countTotT+i)=k;
            tangent_rPO(4,countTotT+i)=gt(k,i);
        end
        clear gt
        countTotT=countTotT+countT;
        
    end
    
    
    
    %Remove the tangents which intersect the obstacles
    countRemoveTang=0;
    indRemoveTang=[];
    for kt=1:countTotT
        flagRemoveTang=0;
        r1=rP;
        r2=tangent_rPO(1:2,kt);
        dr=r2-r1;
        drx=r2(1)-r1(1);
        dry=r2(2)-r1(2);
        mL1=(r2(2)-r1(2))/(r2(1)-r1(1));
        cL1=r1(2)-mL1*r1(1);
        dr_hat=dr/norm(dr);
        %ki1=tangent_rPO{kt}(3,1);
        %ki2=tangent_rPO{kt}(3,2);
        g1=tangent_rPO(4,kt);
        %g2=tangent_rPO{kt}(4,2);
        for k=1:NO
            %if k~=ki1 && k~=ki2
            %countCr=0;
            %check with the circles at the inner vertices
            for i=1:nVO(k)
                P1V=rVO{k}(1:2,i)-r1;
                distLV=abs(dr_hat(1)*P1V(2)-dr_hat(2)*P1V(1));
                if distLV<rho_safe*0.99999999  %if lines intersects the circle at the vertex
                    %tangent_rPO(kt)=[];
                    b=-2*(rVO{k}(1,i)-mL1*cL1+mL1*rVO{k}(2,i));
                    a=(1+mL1^2);
                    c=rVO{k}(1,i)^2+cL1^2+rVO{k}(2,i)^2-2*cL1*rVO{k}(2,i)-rho_safe^2;
                    x_int=(-b+sqrt(b^2-4*a*c))/(2*a);
                    lambda1=(x_int-r1(1))/drx;
                    if (0< lambda1) && (lambda1<1*0.99999999)  % %if line segments intersects the line joining two vertices
                        countRemoveTang=countRemoveTang+1;
                        indRemoveTang(countRemoveTang)=kt;
                        flagRemoveTang=1;
                        break;
                    end
                end
            end
            if flagRemoveTang==1
                break;
            end
            %check with the lines between the inner vertices
            for i=1:nVO(k)
                mL2=(rVO{k}(2,i)-rVO{k}(2,i+1))/(rVO{k}(1,i)-rVO{k}(1,i+1));
                cL2=rVO{k}(2,i)-mL2*rVO{k}(1,i);
                drx2=rVO{k}(1,i+1)-rVO{k}(1,i);
                x_int=(cL2-cL1)/(mL1-mL2);
                %y_int=mL1*x_int+cL1;
                lambda1=(x_int-r1(1))/drx;
                lambda2=(x_int-rVO{k}(1,i))/drx2;
                if (0<= lambda1) && (lambda1<=1) && (0<= lambda2) && (lambda2<=1) % %if line segments intersects the line joining two vertices
                    countRemoveTang=countRemoveTang+1;
                    indRemoveTang(countRemoveTang)=kt;
                    flagRemoveTang=1;
                    break;
                end
            end
            %end
            if flagRemoveTang==1
                break;
            end
        end
    end
    tangent_rPO(:,indRemoveTang)=[];
    
    %plot the tangents which do not intersect obstacles
    if(flag)
        if flag~=0
            fig=flag;
        else
            fig=1;
        end
        figure(fig)
        hold on;
        for kt=1:size(tangent_rPO,2)
            plot([rP(1),tangent_rPO(1,kt)],[rP(2),tangent_rPO(2,kt)],'--','color',colors{6+flag+rp});
            % text(tangent_rPO(1,kt),tangent_rPO(2,kt),[num2str(kt)],'color',colors{8+rp});
        end
    end
    %Add these new nodes to the graph
    indP=indIF(rp);
    indP2=indIF(2)+NTVP;
    %indIF(rp)=indP;
    NTVP=size(tangent_rPO,2);
    rTO_all(:,indP)=rP;
    for i=1:NTVP
        rTO_all(:,indP2+i)=tangent_rPO(1:2,i);
        d=norm(rP-tangent_rPO(1:2,i));
        G(indP,indP2+i)=d;
        G(indP2+i,indP)=d;
        G_pathType(indP,indP2+i)=3;
        G_pathType(indP2+i,indP)=3;
    end
    
    %Now add the segments on the boundary created due to these new nodes
    for i=1:NTVP
        g0=tangent_rPO(4,i);  %Do not store this yet in the global array
        k=tangent_rPO(3,i);
        indPg=find(obsId_all==k);   %indices of the nodes on k^th object
        indg1=find(gTO_all(indPg)>=g0);
        %Now store in the global array
        %         gTO_all(indP2+i)=g0;
        %         obsId_all(indP2+i)=k;
        if isempty(indg1)  %connect to the first and last tangent node on the corresponding obstacle
            ig1=indPg(end);
            ig2=indPg(1);
        elseif indg1(1)==1  %connect to the first and last tangent node on the corresponding obstacle
            ig1=indPg(end);
            ig2=indPg(1);
        else %connect to (indg1-1)^th and indg1^th nodes
            ig1=indPg(indg1(1))-1;
            ig2=indPg(indg1(1));
        end
        %distance with first vertex
        if gTO_all(ig1)<g0
            g01=gTO_all(ig1);
            g02=g0;
            G_pathType(ig1,indP2+i)=1;
            G_pathType(indP2+i,ig1)=2;
        else
            g01=g0;
            g02=gTO_all(ig1);
            G_pathType(ig1,indP2+i)=2;
            G_pathType(indP2+i,ig1)=1;
        end
        indG=find(g01<GammaOk{k} & GammaOk{k}<g02);
        if ~isempty(indG)
            d01=PeriSOk{k}(indG(1)-1)*(GammaOk{k}(indG(1))-g01)/(GammaOk{k}(indG(1))-GammaOk{k}(indG(1)-1));
            d02=PeriSOk{k}(indG(end))*(g02-GammaOk{k}(indG(end)))/(GammaOk{k}(indG(end)+1)-GammaOk{k}(indG(end)));
            % d=d1+d2;
            d1=d01+sum(PeriSOk{k}(indG))+d02;
            if gTO_all(ig1)>g0
                d1=PeriO(k)-d1;
                g02=gTO_all(ig1);
                %switch pathtypes
                G_pathType(ig1,indP2+i)=3-G_pathType(ig1,indP2+i);
                G_pathType(indP2+i,ig1)=3-G_pathType(indP2+i,ig1);
            end
        else
            indG1=find(GammaOk{k}>g01);
            d1= PeriSOk{k}(indG1(1)-1)*(g02-g01)/(GammaOk{k}(indG1(1))-GammaOk{k}(indG1(1)-1));
        end
        %distance with second vertex
        if gTO_all(ig2)<g0
            g01=gTO_all(ig2);
            g02=g0;
            G_pathType(ig2,indP2+i)=1;
            G_pathType(indP2+i,ig2)=2;
        else
            g01=g0;
            g02=gTO_all(ig2);
            G_pathType(ig2,indP2+i)=2;
            G_pathType(indP2+i,ig2)=1;
        end
        indG=find(g01<GammaOk{k} & GammaOk{k}<g02);
        if ~isempty(indG) %endpoints belong to different
            d01=PeriSOk{k}(indG(1)-1)*(GammaOk{k}(indG(1))-g01)/(GammaOk{k}(indG(1))-GammaOk{k}(indG(1)-1));
            d02=PeriSOk{k}(indG(end))*(g02-GammaOk{k}(indG(end)))/(GammaOk{k}(indG(end)+1)-GammaOk{k}(indG(end)));
            % d=d1+d2;
            d2=d01+sum(PeriSOk{k}(indG))+d02;
            
            if gTO_all(ig2)<g0
                d2=PeriO(k)-d2;
                G_pathType(ig2,indP2+i)=3-G_pathType(ig2,indP2+i);
                G_pathType(indP2+i,ig2)=3-G_pathType(indP2+i,ig2);
            end
        else
            indG1=find(GammaOk{k}>g01);
            d2= PeriSOk{k}(indG1(1)-1)*(g02-g01)/(GammaOk{k}(indG1(1))-GammaOk{k}(indG1(1)-1));
        end
        G(ig1,indP2+i)=d1;
        G(indP2+i,ig1)=d1;
        
        G(ig2,indP2+i)=d2;
        G(indP2+i,ig2)=d2;
        
        ind1=find(g0<=GammaOk{k});
        obsVertId_all(indP2+i)=ind1(1)/2; %Id of the vertex of the corresponding obstacle corresponding to the given vertex on the tangent graph
        obsVertPos_all(:,indP2+i)=rVO{k}(:,obsVertId_all(indP2+i)); %Position of the vertex of the corresponding obstacle corresponding to the given vertex on the tangent graph
        %
    end
    
    %Now store in the global array
    for i=1:NTVP
        gTO_all(indP2+i)=tangent_rPO(4,i);
        obsId_all(indP2+i)=tangent_rPO(3,i); %Obstalce Id
    end
    tangent_rPO2{rp}=tangent_rPO;
    clear tangent_rPO;
end

%%
%Find the edges between the new nodes due to rI and rF
tangent_rIO=tangent_rPO2{1};
tangent_rFO=tangent_rPO2{2};
NTVI=size(tangent_rIO,2);
NTVF=size(tangent_rFO,2);
for i=1:NTVI
    k=tangent_rIO(3,i);
    g0=tangent_rIO(4,i);
    ind1=find(tangent_rFO(3,:)==k);
    if ~isempty(ind1)
        for j=1:length(ind1)
            ig1=indIF(2)+NTVI+ind1(j);
            if gTO_all(ig1)<g0
                g01=gTO_all(ig1);
                g02=g0;
                G_pathType(indIF(2)+i,indIF(2)+NTVI+ind1(j))=2;
                G_pathType(indIF(2)+NTVI+ind1(j),indIF(2)+i)=1;
            else
                g01=g0;
                g02=gTO_all(ig1);
                G_pathType(indIF(2)+i,indIF(2)+NTVI+ind1(j))=1;
                G_pathType(indIF(2)+NTVI+ind1(j),indIF(2)+i)=2;
            end
            indG=find(g01<GammaOk{k} & GammaOk{k}<g02);
            if ~isempty(indG)
                d01=PeriSOk{k}(indG(1)-1)*(GammaOk{k}(indG(1))-g01)/(GammaOk{k}(indG(1))-GammaOk{k}(indG(1)-1));
                d02=PeriSOk{k}(indG(end))*(g02-GammaOk{k}(indG(end)))/(GammaOk{k}(indG(end)+1)-GammaOk{k}(indG(end)));
                % d=d1+d2;
                d1=d01+sum(PeriSOk{k}(indG))+d02;
            else
                indG1=find(GammaOk{k}>g01);
                d1= PeriSOk{k}(indG1(1)-1)*(g02-g01)/(GammaOk{k}(indG1(1))-GammaOk{k}(indG1(1)-1));
            end
            
            if d1>PeriO(k)/2  %connect via shortest segment along the boundary
                d1=PeriO(k)-d1;
                %switch the pathtypes
                G_pathType(indIF(2)+i,indIF(2)+NTVI+ind1(j))=3-G_pathType(indIF(2)+i,indIF(2)+NTVI+ind1(j));
                G_pathType(indIF(2)+NTVI+ind1(j),indIF(2)+i)=3-G_pathType(indIF(2)+NTVI+ind1(j),indIF(2)+i);
            end
            
            G(indIF(2)+i,indIF(2)+NTVI+ind1(j))=d1;
            G(indIF(2)+NTVI+ind1(j),indIF(2)+i)=d1;
        end
    end
end


if(flag)
    for i=1:size(rTO_all,2)
        text(rTO_all(1,i),rTO_all(2,i),[num2str(i)],'color',colors{6+flag+rp})
    end
end
%Find the shortest path between rI and rF

[pathLength,path] = dijkstra(G,G,NTV0+1,NTV0+2);

%generate the path between the nodes on the shortest path
r_path=[];
Np=30;
countV0=1;

rV_path(:,1)=rTO_all(:,path(1));
for i=1:length(path)-1
    v1=path(i);
    v2=path(i+1);
    r_path0=[];
    if G_pathType(v1,v2)==3  %Tangent segments
        r_path0=[rTO_all(:,v1),rTO_all(:,v2)];
        k=obsId_all(v2);
    else  %Boundary segments
        g1=gTO_all(v1);
        g2=gTO_all(v2);
        k=obsId_all(v2);
        if G_pathType(v1,v2)==1       %1 for anticlockwise segment on the boundary,
            if g1<=g2
                gamma=g1:(g2-g1)/Np:g2;
                ind1=find(g1<GammaOk{k} & GammaOk{k}<g2);
            else
                gamma=g2:(1-g2)*2/Np:1;
                gamma=[gamma,0:2*g1/Np:g1];
                
                ind1=[find(g2<GammaOk{k} & GammaOk{k}<gammaMax) find(0<=GammaOk{k} & GammaOk{k}<g1)];
            end
        elseif G_pathType(v1,v2)==2  %2 for clockwise segment on the boundary,
            if g1>=g2
                gamma=g1:(g2-g1)/Np:g2;
                ind1=find(g2<GammaOk{k} & GammaOk{k}<g1);
            else
                gamma=g1:-2*g1/Np:0;
                gamma=[gamma,1:(g2-1)*2/Np:g2];
                ind1=[find(g2<GammaOk{k} & GammaOk{k}<gammaMax); find(0<=GammaOk{k} & GammaOk{k}<g1)];
            end
            ind1=flipud(ind1);
        end
        if ~isempty(ind1)
            li=length(ind1);
            rV_path(:,countV0+1:countV0+li)=rVO2{k}(:,ind1);
            rV_path_type(:,countV0+1:countV0+li)=1;  %intersection of straight line and circular ars, (outer vertices)
            rVC_path(:,countV0+1:countV0+li)=rVO{k}(:,ceil(ind1/2));
            pathSegType(1,countV0:countV0+li-1)=ones(1,li)*G_pathType(v1,v2);
            rV_obsId(countV0+1:countV0+li)=k;
        end
        countV0=countV0+length(ind1);
        
        g0=zeros(NO,1);
        for j=1:length(gamma)
            g0(k)=gamma(j);
            r_path0(:,j)=[xOk{k}(g0); yOk{k}(g0)];
        end
    end
    %Add the second vertex on the segment
    pathSegType(countV0)=G_pathType(v1,v2);  %update before counter is updated
    
    countV0=countV0+1;
    rV_path(:,countV0)=rTO_all(:,v2);
    rV_path_type(:,countV0)=0; %Tangent points
    
    rV_obsId(countV0)=k;
    rVC_path(:,countV0)=obsVertPos_all(:,v2);
    r_path=[r_path,r_path0];
end
NrV_path=length(rV_path(1,:));
rV_obsId(NrV_path)=0;  %replace the obsId of the last vertex (final point) by 0
%Remove the intermediate vertices which are redundent  (only keep first, last and intermediate outer vertices on each obstacle in the path)
obsind=find(rV_obsId==rV_obsId(2));
rV_path2=rV_path(:,1);
rVC_path2=rVC_path(:,1);
rV_path_type2=rV_path_type(1);
pathSegType2=pathSegType(1);
while obsind(1)~=1 && obsind(1)~=NrV_path
    ind1=[obsind(1) obsind(find(rV_path_type(obsind)==1)) obsind(end)];
    rV_path2=[rV_path2,rV_path(:,ind1)];
    rVC_path2=[rVC_path2,rVC_path(:,ind1)];
    rV_path_type2=[rV_path_type2,rV_path_type(:,ind1)];
    pathSegType2=[pathSegType2,pathSegType(ind1)];
    %find the indices for the next obstacle
    obsind=find(rV_obsId==rV_obsId(obsind(end)+1));
end
rV_path2=[rV_path2 rV_path(:,end)];
rVC_path2=[rVC_path2 rVC_path(:,end)];
% rV_path_type2=[rV_path_type2 rV_path_type(:,end)];
% pathSegType2=[pathSegType2,pathSegType(end);]
NrV_path2=length(rV_path2(1,:));

%Find the pathlengths
%straight line segments
NVP=length(rV_path2(1,:));
for i=1:2:NVP
    pathLength(i)= norm(rV_path2(:,i)-rV_path2(:,i+1));
end
%cicular segments
for i=2:2:NVP-1
    vec1=rV_path2(:,i)-rVC_path2(:,i);
    vec2=rV_path2(:,i+1)-rVC_path2(:,i);
    deltaAng=acos(vec1'*vec2/norm(vec1)/norm(vec2));
    pathLength(i)= rho_safe*deltaAng;
end

S(1)=0;
for i=2:NVP
    S(i)=S(i-1)+pathLength(i-1);
end

%Return the path and new tangent graph
Path.nodes=path;
Path.rV=rV_path2;
Path.P=pathLength;
Path.S=S;
Path.rVC=rVC_path2;
Path.segType=pathSegType2;
Path.NS=length(pathSegType2);  %Number of segments on the path

tanG_prime.G=G;
tanG_prime.G_pathType=G_pathType;
tanG_prime.rVO=rVO;
tanG_prime.rVO2=rVO2;
tanG_prime.gTO_all=gTO_all;
tanG_prime.rTO_all=rTO_all;
tanG_prime.obsId_all=obsId_all;
tanG_prime.obsVertId_all= obsVertId_all;
tanG_prime.obsVertPos_all= obsVertPos_all;
end
