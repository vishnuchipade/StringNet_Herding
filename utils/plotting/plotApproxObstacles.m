% ------------------------------
% Author: Vishnu S. Chipade
% Email: vishnuc@umich.edu
% ------------------------------

% SetPlotDefaults;
colors={[0,0,1],[1,0,0],[0,0,0],[0.5,0.5,0.5],[0,0.5,1],[1,0.5,1],[0.5,0.5,0],[1,1,0]};
fontSize=18;

% Convex polygon and its centroid
%rVO={[0,5;5,5;4,10;0,12]',[10,-10;15,-10;14,-5]',[-10,-10;-10,-6;-15,-4;-18,-7;-13,-10]'};
rVO={[-10,-9;-9.5,-6;-15,-4;-18,-7;-13,-10]'};
NO=size(rVO,2);
rho_safe=1.5;

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

%Plot the circle around a vertex
ang20=ang2;
ang1=0;
ang2=2*pi;
plot(xV+rho_safe*cos(ang1:(ang2-ang1)/50:ang2),yV+rho_safe*sin(ang1:(ang2-ang1)/50:ang2),'--','color',colors{4})
plot([xV,xV+rho_safe*cos(ang20)],[yV,yV+rho_safe*sin(ang20)],'--','color',colors{4})
text(xV+rho_safe/2*cos(ang20)+0.1,yV+rho_safe/2*sin(ang20),'$\rho_{oa}$','fontsize',fontSize)
