function tanG=tangentGraph(rVO)
% ------------------------------
% Author: Vishnu S. Chipade
% Email: vishnuc@umich.edu
% ------------------------------


%rVO={[0,1;6,1;6,5;0,5;]',[14,3;20,3;20,9;14,9]',[4,-14;11,-14;11,-8;4,-8;]',[7,16;14,16;14,22;7,22;]'};

%rVO={[10,3;14,2;15,7;10.5,8]',[2,-10;7,-12;9,-5;6,-4;]',[12,15;14,18;11,19]'};

global rho_safe 

%rho_safe=1;
options=optimoptions('fmincon','display','off','algorithm','interior-point');%sqp-legacy');
colors={[0,0,1],[1,0,0],[0,0,0],[0.5,0.5,0.5],[0,0.5,1],[1,0.5,1],[0.5,0.5,1],[1,.1,0.5],[0.1,0.9,0.4],[0.9,0.4,0.5],[0,0.5,0.4]};
fontSize=18;

NO=length(rVO);
g=sym('g',[NO,1]);

figure(1)
hold all
%Find the vertices of approximated obstacles

for k=1:NO
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
        
        if(1)
        %plot the circular arc
        plot(xV+rho_safe*cos(ang1:(ang2-ang1)/50:ang2),yV+rho_safe*sin(ang1:(ang2-ang1)/50:ang2),'color',colors{4})
        %plot the straight line
        plot(posVO2(1,2*ii:2*ii+1),posVO2(2,2*ii:2*ii+1),'color',colors{4})
        end
        
    end
    %Find the total perimeter of the boundary of the approximated obstacle
    PeriO(k)=sum(PeriSOk{k});
    GammaOk{k}(1,1)=0;
    for ii=1:nVO(k)
        GammaOk{k}(2*ii,1)=GammaOk{k}(2*ii-1,1)+PeriSOk{k}(2*ii-1)/PeriO(k);
        GammaOk{k}(2*ii+1,1)=GammaOk{k}(2*ii,1)+PeriSOk{k}(2*ii)/PeriO(k);
        GammaCOk{k}(2*ii-1,1)=AngOk{k}(2*ii-1,1)/2/pi;
        GammaCOk{k}(2*ii,1)=AngOk{k}(2*ii,1)/2/pi;
    end
    
    %gk=G(k);
    % syms xOk(gk) yOk(gk)
    xOk{k}=@(g) 0;
    yOk{k}=@(g) 0;
    xCOk{k}=@(g) 0;
    yCOk{k}=@(g) 0;
    dxOk_dg{k}=@(g) 0;
    dyOk_dg{k}=@(g) 0;
    for ii=1:nVO(k)
        i1=2*ii-1;
        i2=2*ii;
        i3=2*ii+1;
        %Get the inner vertex coordinates
        xV=posVO(1,ii);
        yV=posVO(2,ii);
        gk1=GammaOk{k}(i1);
        gk2=GammaOk{k}(i2);
        gk3=GammaOk{k}(i3);
        ang1=AngOk{k}(2*ii-1,1);
        ang2=AngOk{k}(2*ii,1);
        
        if ii==nVO(k)
            xOk{k}=@(g) xOk{k}(g)+(gk1 <= g(k) & g(k) <= gk2)*( xV + rho_safe*cos(ang1+(g(k)-gk1)*(ang2-ang1)/(gk2-gk1)))+...
                (gk2 < g(k) & g(k) <= gk3)*( posVO2(1,i2)+(g(k)-gk2)*(posVO2(1,i3)-posVO2(1,i2))/(gk3-gk2));
            
            yOk{k}=@(g) yOk{k}(g)+(gk1 <= g(k) & g(k) <= gk2)*( yV + rho_safe*sin(ang1+(g(k)-gk1)*(ang2-ang1)/(gk2-gk1)))+...
                (gk2 < g(k) & g(k) <= gk3)*( posVO2(2,i2)+(g(k)-gk2)*(posVO2(2,i3)-posVO2(2,i2))/(gk3-gk2));
            
            xCOk{k}=@(g) xCOk{k}(g)+( GammaCOk{k}(i1) <= g(k) & g(k) <= GammaCOk{k}(i2))*( xV + rho_safe*cos(2*pi*g(k)));
            
            yCOk{k}=@(g) yCOk{k}(g)+(GammaCOk{k}(i1) <= g(k) & g(k) <= GammaCOk{k}(i2))*( yV + rho_safe*sin(2*pi*g(k)));
            
            dxOk_dg{k}=@(g) dxOk_dg{k}(g)+(gk1 <= g(k) & g(k) <= gk2)*(-rho_safe*sin(ang1+(g(k)-gk1)*(ang2-ang1)/(gk2-gk1)))*(ang2-ang1)/(gk2-gk1)+...
                (gk2 < g(k) & g(k) <= gk3)*((posVO2(1,i3)-posVO2(1,i2))/(gk3-gk2));
            
            dyOk_dg{k}=@(g) dyOk_dg{k}(g)+(gk1 <= g(k) & g(k) <= gk2)*(rho_safe*cos(ang1+(g(k)-gk1)*(ang2-ang1)/(gk2-gk1)))*(ang2-ang1)/(gk2-gk1)+...
                (gk2 < g(k) & g(k) <= gk3)*((posVO2(2,i3)-posVO2(2,i2))/(gk3-gk2));
        else
            xOk{k}=@(g) xOk{k}(g)+(gk1 <= g(k) & g(k) <= gk2)*( xV + rho_safe*cos(ang1+(g(k)-gk1)*(ang2-ang1)/(gk2-gk1)))+...
                (gk2 < g(k) & g(k) < gk3)*( posVO2(1,i2)+(g(k)-gk2)*(posVO2(1,i3)-posVO2(1,i2))/(gk3-gk2));
            
            yOk{k}=@(g) yOk{k}(g)+(gk1 <= g(k) & g(k) <= gk2)*( yV + rho_safe*sin(ang1+(g(k)-gk1)*(ang2-ang1)/(gk2-gk1)))+...
                (gk2 < g(k) & g(k) < gk3)*( posVO2(2,i2)+(g(k)-gk2)*(posVO2(2,i3)-posVO2(2,i2))/(gk3-gk2));
            
            xCOk{k}=@(g) xCOk{k}(g)+( GammaCOk{k}(i1) <= g(k) & g(k) < GammaCOk{k}(i2))*( xV + rho_safe*cos(2*pi*g(k)));
            
            yCOk{k}=@(g) yCOk{k}(g)+(GammaCOk{k}(i1) <= g(k) & g(k) < GammaCOk{k}(i2))*( yV + rho_safe*sin(2*pi*g(k)));
            
            dxOk_dg{k}=@(g) dxOk_dg{k}(g)+(gk1 <= g(k) & g(k) <= gk2)*(-rho_safe*sin(ang1+(g(k)-gk1)*(ang2-ang1)/(gk2-gk1)))*(ang2-ang1)/(gk2-gk1)+...
                (gk2 < g(k) & g(k) < gk3)*((posVO2(1,i3)-posVO2(1,i2))/(gk3-gk2));
            
            dyOk_dg{k}=@(g) dyOk_dg{k}(g)+(gk1 <= g(k) & g(k) <= gk2)*(rho_safe*cos(ang1+(g(k)-gk1)*(ang2-ang1)/(gk2-gk1)))*(ang2-ang1)/(gk2-gk1)+...
                (gk2 < g(k) & g(k) < gk3)*((posVO2(2,i3)-posVO2(2,i2))/(gk3-gk2));
        end
        
        Px{k}(ii)=(gk1 <= g(k) & g(k) < gk2)*( xV + rho_safe*cos(ang1+(g(k)-gk1)*(ang2-ang1)/(gk2-gk1)))+...
            (gk2 <= g(k) & g(k) < gk3)*( posVO2(1,i2)+(g(k)-gk2)*(posVO2(1,i3)-posVO2(1,i2))/(gk3-gk2));
        dPx_dg{k}(ii)=(gk1 <= g(k) & g(k) < gk2)*(-rho_safe*sin(ang1+(g(k)-gk1)*(ang2-ang1)/(gk2-gk1)))*(ang2-ang1)/(gk2-gk1)+...
            (gk2 <= g(k) & g(k) < gk3)*((posVO2(1,i3)-posVO2(1,i2))/(gk3-gk2));
        
        d2Px_dg2{k}(ii)=(gk1 <= g(k) & g(k) < gk2)*(-rho_safe*cos(ang1+(g(k)-gk1)*(ang2-ang1)/(gk2-gk1)))*((ang2-ang1)/(gk2-gk1))^2;
        
        Py{k}(ii)=(gk1 <= g(k) & g(k) < gk2)*( yV + rho_safe*sin(ang1+(g(k)-gk1)*(ang2-ang1)/(gk2-gk1)))+...
            (gk2 <= g(k) & g(k) < gk3)*( posVO2(2,i2)+(g(k)-gk2)*(posVO2(2,i3)-posVO2(2,i2))/(gk3-gk2));
        
        dPy_dg{k}(ii)=(gk1 <= g(k) & g(k) < gk2)*(rho_safe*cos(ang1+(g(k)-gk1)*(ang2-ang1)/(gk2-gk1)))*(ang2-ang1)/(gk2-gk1)+...
            (gk2 <= g(k) & g(k) < gk3)*((posVO2(2,i3)-posVO2(2,i2))/(gk3-gk2));
        %xOk=matlabfunction(xOk,'Vars',gk);
        d2Py_dg2{k}(ii)=(gk1 <= g(k) & g(k) < gk2)*(-rho_safe*sin(ang1+(g(k)-gk1)*(ang2-ang1)/(gk2-gk1)))*((ang2-ang1)/(gk2-gk1))^2;
    end
    
    clear posVO2;
    clear posVO;
end
%%
%Find the tangents betweeen all the pairs of the obstacles
A=[eye(NO);-eye(NO)];
%b=[1,1,0,0]';

gamma=0:.01:1;
countP=0;
countTotT=0;
%countTO=zeros(NO,1);
for k=1:NO-1
    for kk=k+1:NO
        countP=countP+1;
        m=@(g) (yOk{k}(g) - yOk{kk}(g))/(xOk{k}(g) - xOk{kk}(g));
        m1=@(g) dyOk_dg{k}(g)/dxOk_dg{k}(g);
        m2=@(g) dyOk_dg{kk}(g)/dxOk_dg{kk}(g);
        f=@(g) (m(g)-m1(g))^2+(m(g)-m2(g))^2;
        
        mC=@(g)  (yCOk{k}(g) - yCOk{kk}(g))/(xCOk{k}(g) - xCOk{kk}(g));
        mC1=@(g) tan(2*pi*g(k)+pi/2);
        mC2=@(g) tan(2*pi*g(kk)+pi/2);
        fC=@(g) (mC(g)-mC1(g))^2+(mC(g)-mC2(g))^2;
        fC=@(g) (tanh(mC(g))-tanh(mC1(g)))^2+(tanh(mC(g))-tanh(mC2(g)))^2;
        
        countT=0;
        
        for i=1:nVO(k)
            for j=1:nVO(kk)
                b=zeros(2*NO,1);
                b([k,kk],1)=[GammaOk{k}(2*i), GammaOk{kk}(2*j)]';
                b(NO+[k,kk],1)=[-GammaOk{k}(2*i-1), -GammaOk{kk}(2*j-1)]';
                g0=zeros(1,NO);
                alp1=.85;
                g0(k)=alp1*GammaOk{k}(2*i-1)+(1-alp1)*GammaOk{k}(2*i);
                g0(kk)=alp1*GammaOk{kk}(2*j-1)+(1-alp1)*GammaOk{kk}(2*j);
                gt0=fmincon(f,g0,A,b,[],[],[],[],[],options);
                ftol=1e-7;
                gtol=1e-5;
                % gt0=fmincon(f,[GammaOk{k}(2*i-1),GammaOk{kk}(2*j-1)],A,b);
                if f(gt0)<ftol && GammaOk{k}(2*i-1)<gt0(k) && gt0(k)<GammaOk{k}(2*i) && GammaOk{kk}(2*j-1)<gt0(kk) && gt0(kk)<GammaOk{kk}(2*j)
                    if countT==0
                        countT=countT+1;
                        gt(:,countT)=gt0;
                    elseif ~sum( (abs(gt(k,:)-gt0(k))<gtol) & (abs(gt(kk,:)-gt0(kk))<gtol))
                        countT=countT+1;
                        gt(:,countT)=gt0;
                    end
                end
                alp2=0.1;
                g0(k)=alp2*GammaOk{k}(2*i-1)+(1-alp2)*GammaOk{k}(2*i);
                g0(kk)=alp2*GammaOk{kk}(2*j-1)+(1-alp2)*GammaOk{kk}(2*j);
                gt0=fmincon(f,g0,A,b,[],[],[],[],[],options);
                % gt0=fmincon(f,[GammaOk{k}(2*i-1),GammaOk{kk}(2*j-1)],A,b);
                if f(gt0)<ftol && GammaOk{k}(2*i-1)<gt0(k) && gt0(k)<GammaOk{k}(2*i) && GammaOk{kk}(2*j-1)<gt0(kk) && gt0(kk)<GammaOk{kk}(2*j)
                    if countT==0
                        countT=countT+1;
                        gt(:,countT)=gt0;
                    elseif ~sum( (abs(gt(k,:)-gt0(k))<gtol) & (abs(gt(kk,:)-gt0(kk))<gtol))
                        countT=countT+1;
                        gt(:,countT)=gt0;
                    end
                end
                
                g0(k)=alp1*GammaOk{k}(2*i-1)+(1-alp1)*GammaOk{k}(2*i);
                g0(kk)=alp2*GammaOk{kk}(2*j-1)+(1-alp2)*GammaOk{kk}(2*j);
                gt0=fmincon(f,g0,A,b,[],[],[],[],[],options);
                % gt0=fmincon(f,[GammaOk{k}(2*i-1),GammaOk{kk}(2*j-1)],A,b);
                if f(gt0)<ftol && GammaOk{k}(2*i-1)<gt0(k) && gt0(k)<GammaOk{k}(2*i) && GammaOk{kk}(2*j-1)<gt0(kk) && gt0(kk)<GammaOk{kk}(2*j)
                    if countT==0
                        countT=countT+1;
                        gt(:,countT)=gt0;
                    elseif ~sum( (abs(gt(k,:)-gt0(k))<gtol) & (abs(gt(kk,:)-gt0(kk))<gtol))
                        countT=countT+1;
                        gt(:,countT)=gt0;
                    end
                end
                
                g0(k)=alp2*GammaOk{k}(2*i-1)+(1-alp2)*GammaOk{k}(2*i);
                g0(kk)=alp1*GammaOk{kk}(2*j-1)+(1-alp1)*GammaOk{kk}(2*j);
                gt0=fmincon(f,g0,A,b,[],[],[],[],[],options);
                % gt0=fmincon(f,[GammaOk{k}(2*i-1),GammaOk{kk}(2*j-1)],A,b);
                if f(gt0)<ftol && GammaOk{k}(2*i-1)<gt0(k) && gt0(k)<GammaOk{k}(2*i) && GammaOk{kk}(2*j-1)<gt0(kk) && gt0(kk)<GammaOk{kk}(2*j)
                    if countT==0
                        countT=countT+1;
                        gt(:,countT)=gt0;
                    elseif ~sum( (abs(gt(k,:)-gt0(k))<gtol) & (abs(gt(kk,:)-gt0(kk))<gtol))
                        countT=countT+1;
                        gt(:,countT)=gt0;
                    end
                end
                
            end
        end
        %X=-2:0.1:15;
        for i=1:countT
            r1(:,i)=[xOk{k}([gt(:,i)']),yOk{k}(gt(:,i)')];
            r2(:,i)=[xOk{kk}([gt(:,i)']),yOk{kk}(gt(:,i)')];
            
            % Store the tangent points on the curves
%             countTO(k)= countTO(k)+1;
%             countTO(kk)= countTO(kk)+1;
%             rTO{k}(:,countTO(k))=r1(:,i);
%             rTO{kk}(:,countTO(kk))=r2(:,i);
%             gTO{k}(countTO(k))=gt(k,i);
%             gTO{kk}(countTO(kk))=gt(kk,i);
            
            %Store the information of tangents globally
            tangentO{countTotT+i}(1:2,1)=r1(:,i);
            tangentO{countTotT+i}(1:2,2)=r2(:,i);
            tangentO{countTotT+i}(3,1:2)=[k,kk];
            tangentO{countTotT+i}(4,1:2)=[gt(k,i),gt(kk,i)];
        end
        hold on;
        for i=1:size(gt,2)
           % plot([r1(1,i),r2(1,i)],[r1(2,i),r2(2,i)],'color',colors{countP+4});
        end
        
        clear gt
        countTotT=countTotT+countT;
    end
    
end


%Remove the tangents which intersect the obstacles
countRemoveTang=0;
indRemoveTang=[];
%tangentO=[];
for kt=1:countTotT
    flagRemoveTang=0;
    r1=tangentO{kt}(1:2,1);
    r2=tangentO{kt}(1:2,2);
    dr=r2-r1;
    drx=r2(1)-r1(1);
    dry=r2(2)-r1(2);
    mL1=(r2(2)-r1(2))/(r2(1)-r1(1));
    cL1=r1(2)-mL1*r1(1);
    dr_hat=dr/norm(dr);
    ki1=tangentO{kt}(3,1);
    ki2=tangentO{kt}(3,2);
    g1=tangentO{kt}(4,1);
    g2=tangentO{kt}(4,2);
    for k=1:NO
        if k~=ki1 && k~=ki2
            %countCr=0;
            %check with the circles at the inner vertices
            for i=1:nVO(k)
                P1V=rVO{k}(1:2,i)-r1;
                distLV=abs(dr_hat(1)*P1V(2)-dr_hat(2)*P1V(1));
                if distLV<rho_safe  %if lines intersects the circle at the vertex
                    %tangentO(kt)=[];
                    b=-2*(rVO{k}(1,i)-mL1*cL1+mL1*rVO{k}(2,i));
                    a=(1+mL1^2);
                    c=rVO{k}(1,i)^2+cL1^2+rVO{k}(2,i)^2-2*cL1*rVO{k}(2,i)-rho_safe^2;
                    x_int=(-b+sqrt(b^2-4*a*c))/(2*a);
                    lambda1=(x_int-r1(1))/drx;
                    if (0<= lambda1) && (lambda1<=1)  % %if line segments intersects the line joining two vertices
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
        end
        if flagRemoveTang==1
            break;
        end
    end
end
tangentO(indRemoveTang)=[];

%%
%Find edge lengths between the nodes on the tangent graph
NTO=length(tangentO);  %Number of tangents
countTO=zeros(NO,1);
for kt=1:NTO    
    r1=tangentO{kt}(1:2,1);
    r2=tangentO{kt}(1:2,2);
    dr=r2-r1;
    drx=r2(1)-r1(1);
    dry=r2(2)-r1(2);
    mL1=(r2(2)-r1(2))/(r2(1)-r1(1));
    cL1=r1(2)-mL1*r1(1);
    dr_hat=dr/norm(dr);
    ki1=tangentO{kt}(3,1);
    ki2=tangentO{kt}(3,2);
    g1=tangentO{kt}(4,1);
    g2=tangentO{kt}(4,2);
    
    % Store the tangent points on the curves
    countTO(ki1)= countTO(ki1)+1;
    countTO(ki2)= countTO(ki2)+1;
    rTO{ki1}(:,countTO(ki1))=r1;
    rTO{ki2}(:,countTO(ki2))=r2;
    gTO{ki1}(countTO(ki1))=g1;
    gTO{ki2}(countTO(ki2))=g2;
    tangentIdO{ki1}(countTO(ki1))=kt;
    tangentIdO{ki2}(countTO(ki2))=kt;
end

gTO_all=[];
rTO_all=[];
tangentIdO_all=[];
NTV=0;  %Total number of tangent vertices
for k=1:NO
        %order the points on the boundary of the obstacles in increasing
        %order of the gamma value (and number them in ascending order)
        [gTO{k} indTO]=sort(gTO{k});
        rTO{k}=rTO{k}(:,indTO);
        tangentIdO{k}=tangentIdO{k}(indTO);
        
        %Find edge lengths along the boundary between the subsequent tangent points
        NTVO=length(gTO{k});
        for kg=1:NTVO-1
            indG=find(gTO{k}(kg)<GammaOk{k} & GammaOk{k}<gTO{k}(kg+1));
            if ~isempty(indG)   %endpoints belong to different circular arcs
            d1=PeriSOk{k}(indG(1)-1)*(GammaOk{k}(indG(1))-gTO{k}(kg))/(GammaOk{k}(indG(1))-GammaOk{k}(indG(1)-1));
            d2=PeriSOk{k}(indG(end))*(gTO{k}(kg+1)-GammaOk{k}(indG(end)))/(GammaOk{k}(indG(end)+1)-GammaOk{k}(indG(end)));
           % d=d1+d2;
            DistTVOk{k}(kg)=d1+sum(PeriSOk{k}(indG(1:end-1)))+d2;
            else  %belong to the same circular arcs
                indG1=find(GammaOk{k}>gTO{k}(kg+1));
               DistTVOk{k}(kg)= PeriSOk{k}(indG1(1)-1)*(gTO{k}(kg+1)-gTO{k}(kg))/(GammaOk{k}(indG1(1))-GammaOk{k}(indG1(1)-1));
            end
            
            G(NTV+kg,NTV+kg+1)= DistTVOk{k}(kg);
            G(NTV+kg+1,NTV+kg)= DistTVOk{k}(kg);
            G_pathType(NTV+kg,NTV+kg+1)= 1;  %1 for anticlockwise segment on the boundary, 2 for a clockwise segment and 3 for tangent segment
            G_pathType(NTV+kg+1,NTV+kg)= 2;
            
            ind1=find(gTO{k}(kg)<=GammaOk{k})
            obsVertId_all(NTV+kg)=ind1(1)/2; %Id of the vertex of the corresponding obstacle corresponding to the given vertex on the tangent graph
            obsVertPos_all(:,NTV+kg)=rVO{k}(:,obsVertId_all(NTV+kg)); %Position of the vertex of the corresponding obstacle corresponding to the given vertex on the tangent graph
%             
        end
        %For the first and the last tangent point
        DistTVOk{k}(NTVO)=PeriO(k)-sum(DistTVOk{k}(1:NTVO-1));
        G(NTV+1,NTV+NTVO)= DistTVOk{k}(NTVO);
        G(NTV+NTVO,NTV+1)= DistTVOk{k}(NTVO);
        G_pathType(NTV+1,NTV+NTVO)= 2;   
        G_pathType(NTV+NTVO,NTV+1)= 1;
        
        ind1=find(gTO{k}(NTVO)<=GammaOk{k})
            obsVertId_all(NTV+NTVO)=ind1(1)/2; %Id of the vertex of the corresponding obstacle corresponding to the given vertex on the tangent graph
            obsVertPos_all(:,NTV+NTVO)=rVO{k}(:,obsVertId_all(NTV+NTVO)); %Position of the vertex of the corresponding obstacle corresponding to the given vertex on the tangent graph
%             
        
        %Collect all the tangents on all the obstacles together
        gTO_all=[gTO_all,gTO{k}];
        rTO_all=[rTO_all,rTO{k}];
        tangentIdO_all=[tangentIdO_all,tangentIdO{k}];
        obsId_all(NTV+1:NTV+NTVO,1)=k; 
        vertId_all(NTV+1:NTV+NTVO,1)=[NTV+1:NTV+NTVO]';  %Id of the vertex of the tangent graph
        
        NTV=NTV+NTVO;
end

[tangentIdO_temp, indT]=sort(tangentIdO_all);
%obsId_temp=obsId_all(indT);
rTO_temp=rTO_all(:,indT);
vertId_temp=vertId_all(indT,1);
for kt=1:NTO
    d=norm(rTO_temp(:,2*kt-1)-rTO_temp(:,2*kt));
    G(vertId_temp(2*kt-1),vertId_temp(2*kt))=d;
    G(vertId_temp(2*kt),vertId_temp(2*kt-1))=d;
    G_pathType(vertId_temp(2*kt-1),vertId_temp(2*kt))=3;
    G_pathType(vertId_temp(2*kt),vertId_temp(2*kt-1))=3;
end

%%
%Plot the tangent graph
if (1)
figure(2)
hold all
%Plot the obstacles
for k=1:NO
    nVO(k)=size(rVO{k},2)-1;
    posVO=rVO{k};
    plot(posVO(1,:),posVO(2,:),'color',colors{3})
    posVO2=rVO2{k};
    posVO2=[posVO2,posVO2(:,1)];
    for ii=1:nVO(k)
        xV=posVO(1,ii);
        yV=posVO(2,ii);
        ang1= AngOk{k}(2*ii-1,1);
        ang2= AngOk{k}(2*ii,1);
        %plot the circular arc
        plot(xV+rho_safe*cos(ang1:(ang2-ang1)/50:ang2),yV+rho_safe*sin(ang1:(ang2-ang1)/50:ang2),'color',colors{4})
        %plot the straight line
        plot(posVO2(1,2*ii:2*ii+1),posVO2(2,2*ii:2*ii+1),'color',colors{4})
    end
end

%plot the tangents which do not intersect obstacles
for kt=1:length(tangentO)
    plot(tangentO{kt}(1,:),tangentO{kt}(2,:),'--','color',colors{1});
end
end
%Return the graph
%G=gTO;

%Return a struct with following properties
tanG.NO=NO;
tanG.nVO=nVO;
tanG.xOk=xOk;
tanG.yOk=yOk;
tanG.dxOk_dg=dxOk_dg;
tanG.dyOk_dg=dyOk_dg;
tanG.GammaOk=GammaOk;
tanG.PeriSOk=PeriSOk;
tanG.PeriO=PeriO;
tanG.G=G;
tanG.G_pathType=G_pathType;
tanG.rVO=rVO;
tanG.rVO2=rVO2;
tanG.gTO_all=gTO_all;
tanG.rTO_all=rTO_all;
tanG.obsId_all=obsId_all;
tanG.vertId_all=vertId_all;
tanG.obsVertId_all= obsVertId_all;
tanG.obsVertPos_all= obsVertPos_all;

end