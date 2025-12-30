%this function plots the beta agent position and velocity on the
% ------------------------------
% Author: Vishnu S. Chipade
% Email: vishnuc@umich.edu
% ------------------------------

%continuously differentiable boundary defined piecewise around the
%obstacles

%options = optimset('Display','off','MaxIter',1000);
options=optimoptions('fmincon','display','off','algorithm','interior-point');%sqp-legacy');
colors={[1,0,0],[0,0,1],[0,0,0],[0.5,0.5,0.5],[0,0.5,1],[1,0.5,1],[0.5,0.5,1],[1,.1,0.5],[0.1,0.9,0.4],[0.9,0.4,0.5],[0,0.5,0.4]};
fontSize=18;
rho_safe=1;

rI=[0,-10]';
%rS=[-3,39]';
%rO=[6,21;8,-7]';
rS=[5,7]';
rO=[0,0;8,-7]';
ROS=norm(rS-rO(:,1));

rotSense=[2,1,2];
RO=2;
w=sqrt(2)*RO;
w=4;
h=w*0.75;
RO_max=6;
k=1;
rA=[14,11]';
vA=[-2,4]';
rA0=[5,-1]';
vA0=[-3,4]';
clear  X Y
xl=10;
[X,Y]=meshgrid(-xl+rO(1,1):.5:xl+rO(1,2),-xl+rO(2,2):.5:xl+rO(2,1));

%rVO={[60,310;160,310;160,390;60,390;]',[-300,350;-190,350;-190,560;-300,560]',[70,-280;120,-280;120,-140;70,-140;]',[-230,-160;-100,-160;-100,-60;-230,-60;]'};
rVO={[0,5;5,5;4,13;0,11]',[10,-8;15,-7;14,-1]',[-5,-9;-5,-4;-10,-2;-13,-5;-10,-8]'};
NO=size(rVO,2);

figure
axis equal
hold all;
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
        angVO2(:,countV)=atan2(posVO2(2,countV)-rCO(2,k),posVO2(1,countV)-rCO(1,k));
        
        if angVO2(:,countV)<0
            angVO2(:,countV)= angVO2(:,countV)+2*pi;
        end
        %corresponding to the next vertex
        countV=countV+1;
        posVO2(:,countV)=[xip,yip]';
        angVO2(:,countV)=atan2(posVO2(2,countV)-rCO(2,k),posVO2(1,countV)-rCO(1,k));
        if angVO2(:,countV)<0
            angVO2(:,countV)= angVO2(:,countV)+2*pi;
        end
    end
    %Shift the outer vertex corresponding to the first inner vertex to
    %appropriate position
    posVO2=posVO2(:,[countV,1:countV-1]);
    angVO2=angVO2(:,[countV,1:countV-1]);
    rVO2{k}=posVO2;
    aVO2{k}=angVO2;
    nV2=countV;
    posVO2=[posVO2,posVO2(:,1)];
    %Find the perimeter and Plot the outer vertices
    PeriO(k)=0;
    for ii=1:nVO(k)
        %Get the inner vertex coordinates
        xV=posVO(1,ii);
        yV=posVO(2,ii);
        if ii~=1
            aVO2{k}(2*ii-1)=aVO2{k}(2*ii-1)-aVO2{k}(1);  %shift with respect the first vertex
            if aVO2{k}(2*ii-1)<0
                aVO2{k}(2*ii-1)=aVO2{k}(2*ii-1)+2*pi;
            end
            aVO2{k}(2*ii)=aVO2{k}(2*ii)-aVO2{k}(1);  %shift with respect the first vertex
            if aVO2{k}(2*ii)<0
                aVO2{k}(2*ii)=aVO2{k}(2*ii)+2*pi;
            end
        else
            aVO2{k}(2*ii)=aVO2{k}(2*ii)-aVO2{k}(1);  %shift with respect the first vertex
            if aVO2{k}(2*ii)<0
                aVO2{k}(2*ii)=aVO2{k}(2*ii)+2*pi;
            end
        end
        
        
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




%Find the projection of rA on obstacle 1
angA=atan2(rA(2)-rCO(2,1),rA(1)-rCO(1,1));
angA=angA-aVO2{1}(1);
if angA<0
    angA=angA+2*pi;
end
indA=find(aVO2{1}<angA);
ind=indA(end);
%the projection of A
if mod(ind,2)==0
    r1=rVO2{1}(:,ind);
    r2=rVO2{1}(:,ind+1);
    dx=r2(1)-r1(1);
    dy=r2(2)-r1(2);
    m=dy/dx;  %slope of line joining the defenders
    c=r1(2)-m*r1(1);
    rAProj(1,1)=(m*rA(2)+rA(1)-m*c)/(1+m^2);
    rAProj(2,1)=m*rAProj(1,1)+c;
    if m<1e16
        lambdaAP=(rAProj(1,1)-r1(1))/dx;
    else
        lambdaAP=(rAProj(2,1)-r1(1))/dy;
    end
    if lambdaAP<0
        rAProj(:,1)=r1;
    elseif lambdaAP>1
        rAProj(:,1)=r1;
    end
    rTP=r1-r2;
    rTP=rTP/norm(rTP);
    if rTP'*vA<0
        rTP=-rTP;
    end
else
    ang=atan2(rA(2)-rVO2{1}(2,ind),rA(1)-rVO2{1}(1,ind));
    rAProj=rVO2{1}(1:2,ind)+rho_safe*[cos(ang);sin(ang)];
    rTP=[-sin(ang);cos(ang)];    
end
%Tangent at the projection point
vAProj=rTP'*vA*rTP;



%position of A
plot(rA(1),rA(2),'o','color',colors{1},'markersize',5);
%Projection
plot(rAProj(1),rAProj(2),'square','color',colors{2});

headSize=2.5;
fontSize=15;
%reference line
plot([rCO(1,1),rCO(1,1)+5],[rCO(2,1),rCO(2)],'b--')
%position of A
plot(rA(1),rA(2),'o','color',colors{1},'markersize',5);
if(0)
    plot(rA0(1),rA0(2),'o','color',colors{1},'markersize',5);
end
%line to rA
plot([rCO(1,1),rA(1)],[rCO(2,1),rA(2)],'color',colors{1})
%velocity of A
quiver(rA(1),rA(2),vA(1),vA(2),0.25,'linewidth',1.5,'color',colors{1},'MaxHeadSize',headSize)
text(rA(1)+0.2,rA(2)+.5,'$\mathbf{v}_a$','fontsize',fontSize)
if(0)
    %initial vA
    quiver(rA0(1),rA0(2),vA0(1),vA0(2),0.25,'linewidth',1.5,'color',colors{1},'MaxHeadSize',headSize)
end
%Projection
plot(rAProj(1),rAProj(2),'square','color',colors{2});
%line to rAProj
plot([rCO(1,1),rAProj(1)],[rCO(2,1),rAProj(2)],'color',colors{2})
%Projected velocity
quiver(rAProj(1),rAProj(2),vAProj(1),vAProj(2),0.25,'linewidth',1.5,'color',colors{2},'MaxHeadSize',headSize)
text(rAProj(1)+0.15,rAProj(2)+.5,'$\mathbf{v}_{\beta}$','fontsize',fontSize)



