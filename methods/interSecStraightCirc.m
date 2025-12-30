function interSec=interSecStraightCirc(r1,r2,mL1,cL1,dr_hat,drx,dry,r21,r22,rVC2,segType2)
% ------------------------------
% Author: Vishnu S. Chipade
% Email: vishnuc@umich.edu
% ------------------------------

%This function gives collision segment for robots moving on a straight
%line segment and a circular arc
global rho_D rho_safe

interSec.Flag=0;

P1V=rVC2(1:2)-r1;
distLV=abs(dr_hat(1)*P1V(2)-dr_hat(2)*P1V(1));
if distLV<=rho_safe + 2*rho_D  %if the robots touch each other somewhere on their paths    
    %find the two intersection points of the inner line and the outer circle
    %choose the inner line
    if mL1<1e16
        if rVC2(2)-mL1*rVC2(1)-cL1>0   %if the center is above the path line
            cL11=cL1+rho_D*sqrt(1+mL1^2);
        else
            cL11=cL1-rho_D*sqrt(1+mL1^2);
        end
        
        b=-2*(rVC2(1)-mL1*cL11+mL1*rVC2(2));
        a=(1+mL1^2);
        c=rVC2(1)^2+cL11^2+rVC2(2)^2-2*cL11*rVC2(2)-(rho_safe+rho_D)^2;
        
        x_int1=(-b-sqrt(b^2-4*a*c))/(2*a);
        x_int2=(-b+sqrt(b^2-4*a*c))/(2*a);
        y_int1=mL1*x_int1+cL11;
        y_int2=mL1*x_int2+cL11;
    else 
        if rVC2(1)-r1(1)>0   %if the center is above the path line
            x_int1=r1(1)+rho_D;
        else
            x_int1=r1(1)-rho_D;
        end
        x_int2=x_int1;
        y_int1=rVC2(2)+sqrt((rho_safe+rho_D)^2-(x_int1-rVC2(1))^2);
        y_int2=rVC2(2)-sqrt((rho_safe+rho_D)^2-(x_int1-rVC2(1))^2);
    end
    ang1=atan2(r21(2)-rVC2(2),r21(1)-rVC2(1));
    if ang1<0
        ang1=ang1+2*pi;
    end
    rotMat=[cos(ang1), sin(ang1);  -sin(ang1), cos(ang1)];
    %rotate all the points such first point becomes x axis
    ri1=rotMat*([x_int1,y_int1]'-rVC2(:));  %Rotate relative to the center of arc
    ri2=rotMat*([x_int2,y_int2]'-rVC2(:));
    r22_prime=rotMat*(r22-rVC2(:));
    angi1=atan2(ri1(2),ri1(1));
    angi2=atan2(ri2(2),ri2(1));
    ang2=atan2(r22_prime(2),r22_prime(1));
    
    %Projection on main path line segment
    x_intp1=(x_int1+mL1*(y_int1-cL1))/(1+mL1^2);
    y_intp1=y_int1+(mL1*x_int1+cL1-y_int1)/(1+mL1^2);
    x_intp2=(x_int2+mL1*(y_int2-cL1))/(1+mL1^2);
    y_intp2=y_int2+(mL1*x_int2+cL1-y_int2)/(1+mL1^2);
    if drx~=0
        lambda1=(x_intp1-r1(1))/drx;
        lambda2=(x_intp2-r1(1))/drx;
    else
        lambda1=(y_intp1-r1(2))/dry;
        lambda2=(y_intp2-r1(2))/dry;
    end
    if lambda1>lambda2  %if lambda1 is greater then the first intersection point is actually the second on the line
        angi0=angi1;
        angi1=angi2;
        angi2=angi0;
        lambda0=lambda1;
        lambda1=lambda2;
        lambda2=lambda0;
    end
    
    if segType2==1
        angcheck1=((0<= angi1) && (angi1<=ang2));
        angcheck2=((0<= angi2) && (angi2<=ang2));
    else
        angcheck1=((ang2<= angi1) && (angi1<=0));
        angcheck2=((ang2<= angi2) && (angi2<=0));
    end
    if ((0<= lambda1) && (lambda1<=1)) || ((0<= lambda2) && (lambda2<=1)) || angcheck1 || angcheck2
        %                 interSec.Pos1(:,1)=[r1(1)+lam11*drx, r1(2)+lam11*dry]'; %first intersection point on path 1
        %                 interSec.Pos1(:,2)=[r1(1)+lam12*drx, r1(2)+lam12*dry]'; %second intersection point on path 1
        %                 interSec.Pbar1=interSec.Pbar1+(lam12-lam11)*L1;
        %                 interSec.S10(1)=Si0+lam11*L1;
        %                 interSec.S10(2)=Si0+lam12*L1;
        flag1=0;
        flag2=0;
        if segType2==1  %anticlockwise segment
            if lambda1<0
                if angi2>ang2 && angcheck1
                    if norm(r1-r22)<2*rho_D
                        flag1=1;
                        rc=r22;
                        ind1=1;
                        ind2=1;
                        rc2=r1;
                    end
                elseif angi1>ang2 && angcheck2
                    if norm(r1-r22)<2*rho_D
                        interSec.Flag=1;   %1 if point intersection
                        interSec.Pos1(:,1)=r1;
                        lam12=min(1,lambda2);
                        interSec.Pos1(:,2)=[r1(1)+lam12*drx, r1(2)+lam12*dry]'; %second intersection point on path 1
                        interSec.Pos2(:,1)=rVC2(:)+(rho_safe)*[cos(angi2+ang1),sin(angi2+ang1)]'; %first intersection point on path 2
                        interSec.Pos2(:,2)=r22; %second intersection point on path 2
                    end
                elseif angi1<0 && angcheck2
                    if norm(r1-r21)<2*rho_D
                        interSec.Flag=1;   %1 if point intersection
                        interSec.Pos1(:,1)=r1;
                        lam12=min(1,lambda2);
                        interSec.Pos1(:,2)=[r1(1)+lam12*drx, r1(2)+lam12*dry]'; %second intersection point on path 1
                        interSec.Pos2(:,2)=rVC2(:)+(rho_safe)*[cos(angi2+ang1),sin(angi2+ang1)]'; %second intersection point on path 2
                        interSec.Pos2(:,1)=r21; %first intersection point on path 2
                    end
                elseif angi2<0 && angcheck1
                    if norm(r1-r21)<2*rho_D
                        flag1=1;
                        rc=r21;
                        ind1=2;
                        ind2=1;
                        rc2=r1;
                    end
                elseif angcheck1 && angcheck2 %both intersection points are on the arc
                    interSec.Flag=1;   %1 if point intersection
                    interSec.Pos1(:,1)=r1;
                    lam12=min(1,lambda2);
                    interSec.Pos1(:,2)=[r1(1)+lam12*drx, r1(2)+lam12*dry]'; %second intersection point on path 1
                    %points on the circular arc;
                    d=2*rho_D+rho_safe;
                    a=(rho_safe^2-(2*rho_D)^2+d^2)/(2*d);
                    h=sqrt(rho_safe^2-a^2);
                    P2=rVC2(:)+a*(rVC2(:)-r1)/d;
                    x_int=P2(1)+h*(r1(2)-rVC2(2))/d;
                    y_int=P2(2)-h*(r1(1)-rVC2(1))/d;
                    ri0=rotMat*([x_int,y_int]'-rVC2(:));
                    angi0=atan2(ri0(2),ri0(1));
                    if ~(0<=angi0 && angi0<ang2)
                        x_int=P2(1)+h*(r1(2)-rVC2(2))/d;
                        y_int=P2(2)-h*(r1(1)-rVC2(1))/d;
                    end
                    if angi1<angi2
                        interSec.Pos2(:,2)=rVC2(:)+(rho_safe)*[cos(angi2+ang1),sin(angi2+ang1)]';
                        interSec.Pos2(:,1)=[x_int,y_int]';
                    else
                        interSec.Pos2(:,1)=rVC2(:)+(rho_safe)*[cos(angi1+ang1),sin(angi1+ang1)]';
                        interSec.Pos2(:,2)=[x_int,y_int]';
                    end
                elseif ~(angcheck1 && angcheck2) %both intersection points are outside the arc
                    %Projection on main path line segment
                    x_intp01=(r21(1)+mL1*(r21(2)-cL1))/(1+mL1^2);
                    y_intp01=r21(2)+(mL1*r21(1)+cL1-r21(2))/(1+mL1^2);
                    x_intp02=(r22(1)+mL1*(r22(2)-cL1))/(1+mL1^2);
                    y_intp02=r22(2)+(mL1*r22(1)+cL1-r22(2))/(1+mL1^2);
                    if drx~=0
                        lambda01=(x_intp01-r1(1))/drx;
                        lambda02=(x_intp02-r1(1))/drx;
                    else
                        lambda01=(y_intp01-r1(2))/dry;
                        lambda02=(y_intp02-r1(2))/dry;
                    end
                    if (0<=lambda01 && lambda01<=1) || (0<=lambda01 && lambda01<=1) || norm(r1-r21)<2*rho_D || norm(r1-r22)<2*rho_D
                        flag2=1;
                    end
                end
            elseif lambda2>1
                if angi1>ang2 && angcheck2
                    if norm(r2-r22)<2*rho_D
                        flag1=1;
                        rc=r22;
                        ind1=1;
                        ind2=2;
                        rc2=r2;
                    end
                elseif angi2>ang2 && angcheck1
                    if norm(r2-r22)<2*rho_D
                        interSec.Flag=1;   %1 if point intersection
                        interSec.Pos1(:,2)=r2;
                        lam11=max(0,lambda1);
                        interSec.Pos1(:,1)=[r1(1)+lam11*drx, r1(2)+lam11*dry]'; %second intersection point on path 1
                        interSec.Pos2(:,1)=rVC2(:)+(rho_safe)*[cos(angi1+ang1),sin(angi1+ang1)]'; %first intersection point on path 2
                        interSec.Pos2(:,2)=r22; %second intersection point on path 2
                    end
                elseif angi2<0 && angcheck1
                    if norm(r2-r21)<2*rho_D
                        interSec.Flag=1;   %1 if point intersection
                        interSec.Pos1(:,2)=r2;
                        lam11=max(0,lambda1);
                        interSec.Pos1(:,1)=[r1(1)+lam11*drx, r1(2)+lam11*dry]'; %second intersection point on path 1
                        interSec.Pos2(:,2)=rVC2(:)+(rho_safe)*[cos(angi1+ang1),sin(angi1+ang1)]'; %second intersection point on path 2
                        interSec.Pos2(:,1)=r21; %first intersection point on path 2
                    end
                elseif angi1<0 && angcheck2
                    if norm(r2-r21)<2*rho_D
                        flag1=1;
                        rc=r21;
                        ind1=2;
                        ind2=2;
                        rc2=r2;
                    end
                elseif angcheck1 && angcheck2 %both intersection points are on the arc
                    interSec.Flag=1;   %1 if point intersection
                    interSec.Pos1(:,2)=r2;
                    lam11=max(0,lambda1);
                    interSec.Pos1(:,1)=[r1(1)+lam11*drx, r1(2)+lam11*dry]'; %second intersection point on path 1
                    %points on the circular arc;
                    d=2*rho_D+rho_safe;
                    a=(rho_safe^2-(2*rho_D)^2+d^2)/(2*d);
                    h=sqrt(rho_safe^2-a^2);
                    P2=rVC2(:)+a*(rVC2(:)-r2)/d;
                    x_int=P2(1)+h*(r2(2)-rVC2(2))/d;
                    y_int=P2(2)-h*(r2(1)-rVC2(1))/d;
                    ri0=rotMat*([x_int,y_int]'-rVC2(:));
                    angi0=atan2(ri0(2),ri0(1));
                    if ~(0<=angi0 && angi0<ang2)
                        x_int=P2(1)+h*(r2(2)-rVC2(2))/d;
                        y_int=P2(2)-h*(r2(1)-rVC2(1))/d;
                    end
                    if angi1<angi2
                        interSec.Pos2(:,2)=rVC2(:)+(rho_safe)*[cos(angi2+ang1),sin(angi2+ang1)]';
                        interSec.Pos2(:,1)=[x_int,y_int]';
                    else
                        interSec.Pos2(:,1)=rVC2(:)+(rho_safe)*[cos(angi1+ang1),sin(angi1+ang1)]';
                        interSec.Pos2(:,2)=[x_int,y_int]';
                    end
                elseif ~(angcheck1 && angcheck2)
                    %Projection on main path line segment
                    x_intp01=(r21(1)+mL1*(r21(2)-cL1))/(1+mL1^2);
                    y_intp01=r21(2)+(mL1*r21(1)+cL1-r21(2))/(1+mL1^2);
                    x_intp02=(r22(1)+mL1*(r22(2)-cL1))/(1+mL1^2);
                    y_intp02=r22(2)+(mL1*r22(1)+cL1-r22(2))/(1+mL1^2);
                    if drx~=0
                        lambda01=(x_intp01-r1(1))/drx;
                        lambda02=(x_intp02-r1(1))/drx;
                    else
                        lambda01=(y_intp01-r1(2))/dry;
                        lambda02=(y_intp02-r1(2))/dry;
                    end
                    
                    if (0<=lambda01 && lambda01<=1) || (0<=lambda01 && lambda01<=1) || norm(r2-r21)<2*rho_D || norm(r2-r22)<2*rho_D
                        flag2=1;
                    end
                end
            elseif (0<=lambda1 && lambda1<=1) && (0<=lambda2 && lambda2<=1)
                interSec.Flag=1;   %1 if point intersection
                interSec.Pos1(:,1)=[r1(1)+lambda1*drx, r1(2)+lambda1*dry]'; %first intersection point on path 1
                interSec.Pos1(:,2)=[r1(1)+lambda2*drx, r1(2)+lambda2*dry]'; %second intersection point on path 1
                angfi1=min(angi2,angi1);  %final path intersection points relative to first point
                angfi2=max(angi1,angi2);
                interSec.Pos2(:,1)=rVC2(:)+(rho_safe)*[cos(angfi1+ang1),sin(angfi1+ang1)]'; %first intersection point on path 2
                interSec.Pos2(:,2)=rVC2(:)+(rho_safe)*[cos(angfi2+ang1),sin(angfi2+ang1)]'; %second intersection point on path 2
            end
            
            if flag1==1
                interSec.Flag=1;   %1 if point intersection
                interSec.Pos1(:,ind2)=rc2;
                if mL1<1e16
                    b=-2*(rc(1)-mL1*cL1+mL1*rc(2));
                    a=(1+mL1^2);
                    c=rc(1)^2+cL1^2+rc(2)^2-2*cL1*rc(2)-(2*rho_D)^2;
                    x_int=(-b-sqrt(b^2-4*a*c))/(2*a);
                    y_int=mL1*x_int+cL1;
                    lambda0=(x_int-r1(1))/drx;
                    if ~(0<=lambda0 && lambda0<=1)
                        x_int=(-b+sqrt(b^2-4*a*c))/(2*a);
                        y_int=mL1*x_int+cL1;
                    end
                    interSec.Pos1(:,3-ind2)=[x_int,y_int]';
                else
                    x_int=rc2(1);
                    y_int=rc(2)+sqrt((2*rho_D)^2-(x_int-rc(1))^2);
                    lambda0=(y_int-r1(2))/dry;
                    if ~(0<=lambda0 && lambda0<=1)
                        y_int=rc(2)-sqrt((2*rho_D)^2-(x_int-rc(1))^2);
                    end
                    interSec.Pos1(:,3-ind2)=[x_int,y_int]';
                end
                               
                %points on the circular arc;
                interSec.Pos2(:,3-ind1)=rc;
                %intersection points between two circles
                d=norm(rVC2-rc2);  %Distance between two centers (P0=rVC2,P1=rc2)
                a=(rho_safe^2-(2*rho_D)^2+d^2)/(2*d);
                h=sqrt(rho_safe^2-a^2);
                P2=rVC2(:)+a*(rc2-rVC2(:))/d;
                x_int=P2(1)+h*(rc2(2)-rVC2(2))/d;
                y_int=P2(2)-h*(rc2(1)-rVC2(1))/d;
                ri0=rotMat*([x_int,y_int]'-rVC2(:));
                angi0=atan2(ri0(2),ri0(1));
                if ~(0<=angi0 && angi0<ang2)
                    x_int=P2(1)-h*(rc2(2)-rVC2(2))/d;
                    y_int=P2(2)+h*(rc2(1)-rVC2(1))/d;
                end
                interSec.Pos2(:,ind1)=[x_int,y_int]';
            end
            
        else  %Clockwise segment
            if lambda1<0
                if angi2<ang2 && angcheck1
                    if norm(r1-r22)<2*rho_D
                        flag1=1;
                        rc=r22;
                        ind1=1;
                        ind2=1;
                        rc2=r1;
                    end
                elseif angi1<ang2 && angcheck2
                    if norm(r1-r22)<2*rho_D
                        interSec.Flag=1;   %1 if point intersection
                        interSec.Pos1(:,1)=r1;
                        lam12=min(1,lambda2);
                        interSec.Pos1(:,2)=[r1(1)+lam12*drx, r1(2)+lam12*dry]'; %second intersection point on path 1
                        interSec.Pos2(:,1)=rVC2(:)+(rho_safe)*[cos(angi2+ang1),sin(angi2+ang1)]'; %first intersection point on path 2
                        interSec.Pos2(:,2)=r22; %second intersection point on path 2
                    end
                elseif angi1>0 && angcheck2
                    if norm(r1-r21)<2*rho_D
                        interSec.Flag=1;   %1 if point intersection
                        interSec.Pos1(:,1)=r1;
                        lam12=min(1,lambda2);
                        interSec.Pos1(:,2)=[r1(1)+lam12*drx, r1(2)+lam12*dry]'; %second intersection point on path 1
                        interSec.Pos2(:,2)=rVC2(:)+(rho_safe)*[cos(angi2+ang1),sin(angi2+ang1)]'; %second intersection point on path 2
                        interSec.Pos2(:,1)=r21; %first intersection point on path 2
                    end
                elseif angi2>0 && angcheck1
                    if norm(r1-r21)<2*rho_D
                        flag1=1;
                        rc=r21;
                        ind1=2;
                        ind2=1;
                        rc2=r1;
                    end
                elseif angcheck1 && angcheck2 %both intersection points are on the arc
                    interSec.Flag=1;   %1 if point intersection
                    interSec.Pos1(:,1)=r1;
                    lam12=min(1,lambda2);
                    interSec.Pos1(:,2)=[r1(1)+lam12*drx, r1(2)+lam12*dry]'; %second intersection point on path 1
                    %points on the circular arc;
                    d=2*rho_D+rho_safe;
                    a=(rho_safe^2-(2*rho_D)^2+d^2)/(2*d);
                    h=sqrt(rho_safe^2-a^2);
                    P2=rVC2(:)+a*(rVC2(:)-r1)/d;
                    x_int=P2(1)+h*(r1(2)-rVC2(2))/d;
                    y_int=P2(2)-h*(r1(1)-rVC2(1))/d;
                    ri0=rotMat*([x_int,y_int]'-rVC2(:));
                    angi0=atan2(ri0(2),ri0(1));
                    if ~(0>=angi0 && angi0>ang2)
                        x_int=P2(1)+h*(r1(2)-rVC2(2))/d;
                        y_int=P2(2)-h*(r1(1)-rVC2(1))/d;
                    end
                    if angi1>=angi2
                        interSec.Pos2(:,2)=rVC2(:)+(rho_safe)*[cos(angi2+ang1),sin(angi2+ang1)]';
                        interSec.Pos2(:,1)=[x_int,y_int]';
                    else
                        interSec.Pos2(:,1)=rVC2(:)+(rho_safe)*[cos(angi1+ang1),sin(angi1+ang1)]';
                        interSec.Pos2(:,2)=[x_int,y_int]';
                    end
                elseif ~(angcheck1 && angcheck2) %both intersection points are outside the arc
                    %Projection on main path line segment
                    x_intp01=(r21(1)+mL1*(r21(2)-cL1))/(1+mL1^2);
                    y_intp01=r21(2)+(mL1*r21(1)+cL1-r21(2))/(1+mL1^2);
                    x_intp02=(r22(1)+mL1*(r22(2)-cL1))/(1+mL1^2);
                    y_intp02=r22(2)+(mL1*r22(1)+cL1-r22(2))/(1+mL1^2);
                    if drx~=0
                        lambda01=(x_intp01-r1(1))/drx;
                        lambda02=(x_intp02-r1(1))/drx;
                    else
                        lambda01=(y_intp01-r1(2))/dry;
                        lambda02=(y_intp02-r1(2))/dry;
                    end
                    if (0<=lambda01 && lambda01<=1) || (0<=lambda01 && lambda01<=1) || norm(r1-r21)<2*rho_D || norm(r1-r22)<2*rho_D
                        flag2=1;
                    end
                end
            elseif lambda2>1
                if angi1<ang2 && angcheck2
                    if norm(r2-r22)<2*rho_D
                        flag1=1;
                        rc=r22;
                        ind1=1;
                        ind2=2;
                        rc2=r2;
                    end
                elseif angi2<ang2 && angcheck1
                    if norm(r2-r22)<2*rho_D
                        interSec.Flag=1;   %1 if point intersection
                        interSec.Pos1(:,2)=r2;
                        lam11=max(0,lambda1);
                        interSec.Pos1(:,1)=[r1(1)+lam11*drx, r1(2)+lam11*dry]'; %second intersection point on path 1
                        interSec.Pos2(:,1)=rVC2(:)+(rho_safe)*[cos(angi1+ang1),sin(angi1+ang1)]'; %first intersection point on path 2
                        interSec.Pos2(:,2)=r22; %second intersection point on path 2
                    end
                elseif angi2>0 && angcheck1
                    if norm(r2-r21)<2*rho_D
                        interSec.Flag=1;   %1 if point intersection
                        interSec.Pos1(:,2)=r2;
                        lam11=max(0,lambda1);
                        interSec.Pos1(:,1)=[r1(1)+lam11*drx, r1(2)+lam11*dry]'; %second intersection point on path 1
                        interSec.Pos2(:,2)=rVC2(:)+(rho_safe)*[cos(angi1+ang1),sin(angi1+ang1)]'; %second intersection point on path 2
                        interSec.Pos2(:,1)=r21; %first intersection point on path 2
                    end
                elseif angi1>0 && angcheck2
                    if norm(r2-r21)<2*rho_D
                        flag1=1;
                        rc=r21;
                        ind1=2;
                        ind2=2;
                        rc2=r2;
                    end
                elseif angcheck1 && angcheck2 %both intersection points are on the arc
                    interSec.Flag=1;   %1 if point intersection
                    interSec.Pos1(:,2)=r2;
                    lam11=max(0,lambda1);
                    interSec.Pos1(:,1)=[r1(1)+lam11*drx, r1(2)+lam11*dry]'; %second intersection point on path 1
                    %points on the circular arc;
                    d=2*rho_D+rho_safe;
                    a=(rho_safe^2-(2*rho_D)^2+d^2)/(2*d);
                    h=sqrt(rho_safe^2-a^2);
                    P2=rVC2(:)+a*(rVC2(:)-r2)/d;
                    x_int=P2(1)+h*(r2(2)-rVC2(2))/d;
                    y_int=P2(2)-h*(r2(1)-rVC2(1))/d;
                    ri0=rotMat*([x_int,y_int]'-rVC2(:));
                    angi0=atan2(ri0(2),ri0(1));
                    if ~(0>=angi0 && angi0>ang2)
                        x_int=P2(1)+h*(r2(2)-rVC2(2))/d;
                        y_int=P2(2)-h*(r2(1)-rVC2(1))/d;
                    end
                    if angi1>angi2
                        interSec.Pos2(:,2)=rVC2(:)+(rho_safe)*[cos(angi2+ang1),sin(angi2+ang1)]';
                        interSec.Pos2(:,1)=[x_int,y_int]';
                    else
                        interSec.Pos2(:,1)=rVC2(:)+(rho_safe)*[cos(angi1+ang1),sin(angi1+ang1)]';
                        interSec.Pos2(:,2)=[x_int,y_int]';
                    end
                elseif ~(angcheck1 && angcheck2)
                    %Projection on main path line segment
                    x_intp01=(r21(1)+mL1*(r21(2)-cL1))/(1+mL1^2);
                    y_intp01=r21(2)+(mL1*r21(1)+cL1-r21(2))/(1+mL1^2);
                    x_intp02=(r22(1)+mL1*(r22(2)-cL1))/(1+mL1^2);
                    y_intp02=r22(2)+(mL1*r22(1)+cL1-r22(2))/(1+mL1^2);
                    if drx~=0
                        lambda01=(x_intp01-r1(1))/drx;
                        lambda02=(x_intp02-r1(1))/drx;
                    else
                        lambda01=(y_intp01-r1(2))/dry;
                        lambda02=(y_intp02-r1(2))/dry;
                    end
                    
                    if (0<=lambda01 && lambda01<=1) || (0<=lambda01 && lambda01<=1) || norm(r2-r21)<2*rho_D || norm(r2-r22)<2*rho_D                        
                        flag2=1;
                    end
                end
            elseif (0<=lambda1 && lambda1<=1) && (0<=lambda2 && lambda2<=1)
                interSec.Flag=1;   %1 if point intersection
                interSec.Pos1(:,1)=[r1(1)+lambda1*drx, r1(2)+lambda1*dry]'; %first intersection point on path 1
                interSec.Pos1(:,2)=[r1(1)+lambda2*drx, r1(2)+lambda2*dry]'; %second intersection point on path 1
                angfi1=min(angi2,angi1);  %final path intersection points relative to first point
                angfi2=max(angi1,angi2);
                interSec.Pos2(:,1)=rVC2(:)+(rho_safe)*[cos(angfi1+ang1),sin(angfi1+ang1)]'; %first intersection point on path 2
                interSec.Pos2(:,2)=rVC2(:)+(rho_safe)*[cos(angfi2+ang1),sin(angfi2+ang1)]'; %second intersection point on path 2
            end
            
            if flag1==1
                interSec.Flag=1;   %1 if point intersection
                interSec.Pos1(:,ind2)=rc2;
                if mL1<1e16
                    b=-2*(rc(1)-mL1*cL1+mL1*rc(2));
                    a=(1+mL1^2);
                    c=rc(1)^2+cL1^2+rc(2)^2-2*cL1*rc(2)-(2*rho_D)^2;
                    x_int=(-b-sqrt(b^2-4*a*c))/(2*a);
                    y_int=mL1*x_int+cL1;
                    lambda0=(x_int-r1(1))/drx;
                    if ~(0<=lambda0 && lambda0<=1)
                        x_int=(-b+sqrt(b^2-4*a*c))/(2*a);
                        y_int=mL1*x_int+cL1;
                    end
                    interSec.Pos1(:,3-ind2)=[x_int,y_int]';
                else
                    x_int=rc2(1);
                    y_int=rc(2)+sqrt((2*rho_D)^2-(x_int-rc(1))^2);
                    lambda0=(y_int-r1(2))/dry;
                    if ~(0<=lambda0 && lambda0<=1)
                        y_int=rc(2)-sqrt((2*rho_D)^2-(x_int-rc(1))^2);
                    end
                    interSec.Pos1(:,3-ind2)=[x_int,y_int]';
                end
                %points on the circular arc;
                interSec.Pos2(:,3-ind1)=rc;
                %intersection points between two circles
                d=norm(rVC2-rc2);  %Distance between two centers (P0=rVC2,P1=rc2)
                a=(rho_safe^2-(2*rho_D)^2+d^2)/(2*d);
                h=sqrt(rho_safe^2-a^2);
                P2=rVC2(:)+a*(rc2-rVC2(:))/d;
                x_int=P2(1)+h*(rc2(2)-rVC2(2))/d;
                y_int=P2(2)-h*(rc2(1)-rVC2(1))/d;
                ri0=rotMat*([x_int,y_int]'-rVC2(:));
                angi0=atan2(ri0(2),ri0(1));
                if ~(0>=angi0 && angi0>ang2)
                    x_int=P2(1)+h*(rc2(2)-rVC2(2))/d;
                    y_int=P2(2)-h*(rc2(1)-rVC2(1))/d;
                end
                interSec.Pos2(:,ind1)=[x_int,y_int]';
            end
        end
        
        
        if flag2==1
            interSec.Flag=1;   %1 if point intersection
            interSec.Pos2(:,1)=r21;
            interSec.Pos2(:,2)=r22;
            for ii=1:2
                rc= interSec.Pos2(:,ii);
                if mL1<1e16
                    b=-2*(rc(1)-mL1*cL1+mL1*rc(2));
                    a=(1+mL1^2);
                    c=rc(1)^2+cL1^2+rc(2)^2-2*cL1*rc(2)-(2*rho_D)^2;
                    x_int0(1,ii)=(-b-sqrt(b^2-4*a*c))/(2*a);
                    y_int0(1,ii)=mL1*x_int0(1,ii)+cL1;
                    lambda0(1,ii)=(x_int0(1,ii)-r1(1))/drx;
                    x_int0(2,ii)=(-b+sqrt(b^2-4*a*c))/(2*a);
                    y_int0(2,ii)=mL1*x_int0(2,ii)+cL1;
                    lambda0(2,ii)=(x_int0(2,ii)-r1(1))/drx;
                else
                    x_int0(1,ii)=r1(1);
                    x_int0(2,ii)=r1(1);
                    y_int0(1,ii)=rc(2)+sqrt((2*rho_D)^2-(x_int0(1,ii)-rc(1))^2);
                    lambda0(1,ii)=(y_int0(1,ii)-r1(2))/dry;
                    y_int0(2,ii)=rc(2)-sqrt((2*rho_D)^2-(x_int0(2,ii)-rc(1))^2);
                    lambda0(2,ii)=(y_int0(2,ii)-r1(2))/dry;
                end
            end
            x_int0=x_int0(:);
            y_int0=y_int0(:);
            lambda0=lambda0(:);
            lam11=max(0,min(lambda0));
            lam12=min(1,max(lambda0));
            interSec.Pos1(:,1)=[r1(1)+lam11*drx, r1(2)+lam11*dry]'; %first intersection point on path 1
            interSec.Pos1(:,2)=[r1(1)+lam12*drx, r1(2)+lam12*dry]'; %second intersection point on path 1
        end
        %                 interSec.Pos2(:,1)=rVC2(:)+(rho_safe)*[cos(angfi1+ang1),sin(angfi1+ang1)]'; %first intersection point on path 2
        %                 interSec.Pos2(:,2)=rVC2(:)+(rho_safe)*[cos(angfi2+ang1),sin(angfi2+ang1)]'; %second intersection point on path 2
        if  interSec.Flag~=0
            ds1=norm(interSec.Pos1(:,1)-interSec.Pos1(:,2));
            interSec.Pbar1=ds1;
            interSec.S10(1)=norm(interSec.Pos1(:,1)-r1);
            interSec.S10(2)=interSec.S10(1)+ds1;
            
            ds2=rho_safe*acos(1-norm(interSec.Pos2(:,1)-interSec.Pos2(:,2))^2/(2*rho_safe^2));
            interSec.Pbar2=ds2;
            interSec.S20(1)=rho_safe*acos(1-norm(interSec.Pos2(:,1)-r21)^2/(2*rho_safe^2));
            interSec.S20(2)=interSec.S20(1)+ds2;
            return;
        end
    end
end
end