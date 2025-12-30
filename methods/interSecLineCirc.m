function interSec=interSecLineCirc(r1,r2,mL1,cL1,L1,dr_hat,drx,dry,r21,r22,rVC2,segType2)
% ------------------------------
% Author: Vishnu S. Chipade
% Email: vishnuc@umich.edu
% ------------------------------

%This function gives collision segment for robots moving on a straight
%line segment and a circular arc
global rho_D rho_safe dthetai
%interchange the points and sense of rotation if the segment is clockwise 
%(i.e. convert into anticlockwise, check intersection and convert the intersection points back to clockwise, done at the end)
if segType2==2 
    r22_temp=r22;
    r21_temp=r21;
    r21=r22_temp;
    r22=r21_temp;
end

interSec.Flag=0;

P1V=rVC2(1:2)-r1;
distLV=abs(dr_hat(1)*P1V(2)-dr_hat(2)*P1V(1));
if distLV<=rho_safe + 2*rho_D && distLV>=rho_safe*0.999999999 %if the robots can touch each other somewhere on their paths
    %find the two intersection points of the line 2*rho_D far from the main
    %line and the circular arc
    if mL1<1e16
        if rVC2(2)-mL1*rVC2(1)-cL1>0   %if the center is above the path line
            cL11=cL1+2*rho_D*sqrt(1+mL1^2);
        else
            cL11=cL1-2*rho_D*sqrt(1+mL1^2);
        end
        
        b=-2*(rVC2(1)-mL1*cL11+mL1*rVC2(2));
        a=(1+mL1^2);
        c=rVC2(1)^2+cL11^2+rVC2(2)^2-2*cL11*rVC2(2)-(rho_safe)^2;
        bsquare4ac=b^2-4*a*c;
        if abs(bsquare4ac)<1e-10
            bsquare4ac=0;
        end
        x_int1=(-b-sqrt(bsquare4ac))/(2*a);
        x_int2=(-b+sqrt(bsquare4ac))/(2*a);
        y_int1=mL1*x_int1+cL11;
        y_int2=mL1*x_int2+cL11;
    else
        if rVC2(1)-r1(1)>0   %if the center is above the path line
            x_int1=r1(1)+2*rho_D;
        else
            x_int1=r1(1)-2*rho_D;
        end
        x_int2=x_int1;
        y_int1=rVC2(2)+sqrt((rho_safe)^2-(x_int1-rVC2(1))^2);
        y_int2=rVC2(2)-sqrt((rho_safe)^2-(x_int1-rVC2(1))^2);
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
    
    %Find the intersection point of normal through the center to the line
    if mL1<1e16
        x_int0=(rVC2(1)+mL1*rVC2(2)-cL1*mL1)/(1+mL1^2);
        y_int0=mL1*x_int0+cL1;
    else
        x_int0=r1(1);
        y_int0=rVC2(2);
    end
    ri20=rotMat*([x_int0,y_int0]'-rVC2(:));
    angi0=atan2(ri20(2),ri20(1));
    %Lambda parameters for these intersection points
    if drx~=0
        lambda0=(x_int0-r1(1))/drx;
        lambda1=(x_intp1-r1(1))/drx;
        lambda2=(x_intp2-r1(1))/drx;
    else
        lambda0=(y_int0-r1(2))/dry;
        lambda1=(y_intp1-r1(2))/dry;
        lambda2=(y_intp2-r1(2))/dry;
    end
    if lambda1>lambda2  %if lambda1 is greater then the first intersection point is actually the second on the line
        angi_temp=angi1;
        angi1=angi2;
        angi2=angi_temp;
        lambda_temp=lambda1;
        lambda1=lambda2;
        lambda2=lambda_temp;
    end
    
    %if segType2==1
    angi1check=((0<= angi1) && (angi1<=ang2));
    angi2check=((0<= angi2) && (angi2<=ang2));
    %     else
    %         angi1check=((ang2<= angi1) && (angi1<=0));
    %         angi2check=((ang2<= angi2) && (angi2<=0));
    %     end
    %   if ((0<= lambda1) && (lambda1<=1)) || ((0<= lambda2) && (lambda2<=1)) || angi1check || angi2check
    %                 interSec.Pos1(:,1)=[r1(1)+lam11*drx, r1(2)+lam11*dry]'; %first intersection point on path 1
    %                 interSec.Pos1(:,2)=[r1(1)+lam12*drx, r1(2)+lam12*dry]'; %second intersection point on path 1
    %                 interSec.Pbar1=interSec.Pbar1+(lam12-lam11)*L1;
    %                 interSec.S10(1)=Si0+lam11*L1;
    %                 interSec.S10(2)=Si0+lam12*L1;
    flag1=0;
    flag2=0;
    %if segType2==1  %anticlockwise segment
    if angi0>=ang2
        if lambda0>0 && lambda0<1            
            rp=r22;
            if drx~=0
                x_p=(rp(1)+mL1*(rp(2)-cL1))/(1+mL1^2);
                lambdap=(x_p-r1(1))/drx;
                drpL1=((cL1+mL1*rp(1)-rp(2))/sqrt(1+mL1^2));
            else
                y_p=rp(2);
                lambdap=(y_p-r1(2))/dry;
                drpL1=abs(r1(1)-rp(1));
            end
            dlambdai=2*rho_D/L1;
            %check1=-dthetai<=angr1 && angr1<=ang2+dthetai;
            check2=(-dlambdai<=lambdap && lambdap<=1+dlambdai);
            
            if abs(cL1+mL1*r22(1)-r22(2))/sqrt(1+mL1^2) <=2*rho_D && check2 %distance between the 2nd end point on the arc and the line
                interSec.Flag=1;
                angfi2=ang2;
                %Find two intersection points of the circle at the
                %end point of arc and segment
                rc=r22;
                lambdac=lambdaInterSecCircLine(rc,2*rho_D,r1,mL1,cL1,drx,dry); %gives lambda values of two intersection points on the line and circle centered at rc (radius 2*rho_D)
                
                %Check which side of the circular arc the point r1 is on
                r1_temp=rotMat*(r1-rVC2);
                angr1=atan2(r1_temp(2),r1_temp(1));
                if angr1<angi0
                    if lambda1<=0
                        lam11=0;
                        lam12=min(1,max(lambdac));
                        %Find two intersection points of the circle at the
                        %end point of segment 2 and segment 1
                        rc=r1;
                        angc=angInterSecCircCirc(rc,2*rho_D,rVC2,rho_safe,rotMat);  %angles of intersection points between two circles with respect to first end point
                        angfi1=max(0,min(angc));
                    else
                        lam12=min(1,max(lambdac));
                        if angi1<=0
                            angfi1=0;
                            lambdac=lambdaInterSecCircLine(r21,2*rho_D,r1,mL1,cL1,drx,dry); %gives lambda values of two intersection points on the line and circle centered at rc (radius 2*rho_D)
                            lam11=max(0,min(lambdac));
                        else
                            lam11=lambda1;
                            angfi1=angi1;
                        end
                    end
                else
                    if lambda2>=1
                        lam12=1;
                        lam11=max(0,min(lambdac));
                        %Find two intersection points of the circle at the
                        %end point of segment 2 and segment 1
                        rc=r2;
                        angc=angInterSecCircCirc(rc,2*rho_D,rVC2,rho_safe,rotMat);  %angles of intersection points between two circles with respect to first end point
                        angfi1=max(0,min(angc));
                    else
                        lam11=max(0,min(lambdac));
                        if angi2<=0
                            angfi1=0;
                            lambdac=lambdaInterSecCircLine(r21,2*rho_D,r1,mL1,cL1,drx,dry); %gives lambda values of two intersection points on the line and circle centered at rc (radius 2*rho_D)
                            lam12=min(1,max(lambdac));
                        else
                            lam12=lambda2;
                            angfi1=angi2;
                        end
                    end
                end
            end
        elseif lambda0<=0
            r1_temp=rotMat*(r1-rVC2);
            angr1=atan2(r1_temp(2),r1_temp(1));
            rp=r22;
            if drx~=0
                x_p=(rp(1)+mL1*(rp(2)-cL1))/(1+mL1^2);
                lambdap=(x_p-r1(1))/drx;
                drpL1=((cL1+mL1*rp(1)-rp(2))/sqrt(1+mL1^2));
            else
                y_p=rp(2);
                lambdap=(y_p-r1(2))/dry;
                drpL1=abs(r1(1)-rp(1));
            end
            dlambdai=2*rho_D/L1;
            check1=-dthetai<=angr1 && angr1<=ang2+dthetai;
            check2=(-dlambdai<=lambdap && lambdap<=1+dlambdai);
            if norm(r1-rVC2)<=rho_safe+2*rho_D && drpL1 <= 2*rho_D  && (check1 || check2)  %norm(r22-r1)<=2*rho_D %|| (0<=angr1 && angr1<=ang2)
                interSec.Flag=1;
                if angr1<angi0
                    angfi2=ang2;
                    lambdac0=lambdaInterSecCircLine(r22,2*rho_D,r1,mL1,cL1,drx,dry); %gives lambda values of two intersection points on the line and circle centered at rc (radius 2*rho_D)
                    lambdac=lambdac0(0<lambdac0 & lambdac0<1);
                    if length(lambdac)==2
                        %lambda2=max(lambdac);
                        lam11=min(lambdac);
                    elseif length(lambdac)==1
                        lam11=0;
                    elseif lambdac0(1)<0 && lambdac0(2)<0
                        angc=angInterSecCircCirc(r1,2*rho_D,rVC2,rho_safe,rotMat);  %angles of intersection points between two circles with respect to first end point
                        angfi2=min(ang2,max(angc));
                        lam11=0;
                    elseif max(lambdac0)>1 && min(lambdac0)<0
                        lam12=1;
                        lam11=0;
                    end
                    
                    if angi2<0
                        angfi1=0;
                        lambdac0=lambdaInterSecCircLine(r21,2*rho_D,r1,mL1,cL1,drx,dry); %gives lambda values of two intersection points on the line and circle centered at rc (radius 2*rho_D)
                        lambdac=lambdac0(0<lambdac0 & lambdac0<1);
                        if length(lambdac)==2
                            %lambda2=max(lambdac);
                            lam12=max(lambdac);
                        elseif length(lambdac)==1
                            lam12=lambdac;
                        elseif lambdac0(1)>1 && lambdac0(2)>1
                            angc=angInterSecCircCirc(r2,2*rho_D,rVC2,rho_safe,rotMat);  %angles of intersection points between two circles with respect to first end point
                            angfi1=max(0,min(angc));
                            lam12=1;
                        elseif max(lambdac0)>1 && min(lambdac0)<0
                            lam12=1;
                            lam11=0;
                        end
                    else
                        angfi1=angi2;
                        lam12=lambda2;
                        if lambda2>1
                            angc=angInterSecCircCirc(r2,2*rho_D,rVC2,rho_safe,rotMat);  %angles of intersection points between two circles with respect to first end point
                            angfi1=max(angi2,min(angc));
                            lam12=1;
                        end
                    end
                else
                    angc=angInterSecCircCirc(r1,2*rho_D,rVC2,rho_safe,rotMat);  %angles of intersection points between two circles with respect to first end point
                    lambdac=lambdaInterSecCircLine(r22,2*rho_D,r1,mL1,cL1,drx,dry); %gives lambda values of two intersection points on the line and circle centered at rc (radius 2*rho_D)
                    lam12=min(1,max(lambdac));
                    angfi1=max(0,min(angc));
                    angfi2=ang2;
                    lam11=0;
                end
            end
        elseif lambda0>=1
            r2_temp=rotMat*(r2-rVC2);
            angr2=atan2(r2_temp(2),r2_temp(1));
            rp=r22;
            if drx~=0
                x_p=(rp(1)+mL1*(rp(2)-cL1))/(1+mL1^2);
                lambdap=(x_p-r1(1))/drx;
                drpL1=((cL1+mL1*rp(1)-rp(2))/sqrt(1+mL1^2));
            else
                y_p=rp(2);
                lambdap=(y_p-r1(2))/dry;
                drpL1=abs(r1(1)-rp(1));
            end
            dlambdai=2*rho_D/L1;
            check1=-dthetai<=angr2 && angr2<=ang2+dthetai;
            check2=(-dlambdai<=lambdap && lambdap<=1+dlambdai);
            if norm(r2-rVC2)<=rho_safe+2*rho_D && drpL1 <= 2*rho_D && (check1 || check2)  %norm(r22-r2)<=2*rho_D %|| (0<=angr2 && angr2<=ang2)
                interSec.Flag=1;
                if angr2<angi0
                    angfi2=ang2;
                    lambdac0=lambdaInterSecCircLine(r22,2*rho_D,r1,mL1,cL1,drx,dry); %gives lambda values of two intersection points on the line and circle centered at rc (radius 2*rho_D)
                    lambdac=lambdac0(0<lambdac0 & lambdac0<1);
                    if length(lambdac)==2
                        %lambda2=max(lambdac);
                        lam12=max(lambdac);
                    elseif length(lambdac)==1
                        lam12=1;
                    elseif lambdac0(1)>1 && lambdac0(2)>1
                        angc=angInterSecCircCirc(r2,2*rho_D,rVC2,rho_safe,rotMat);  %angles of intersection points between two circles with respect to first end point
                        angfi2=min(ang2,max(angc));
                        lam12=1;
                    elseif max(lambdac0)>1 && min(lambdac0)<0
                        lam12=1;
                        lam11=0;
                    end
                    
                    if angi1<0
                        angfi1=0;
                        lambdac0=lambdaInterSecCircLine(r21,2*rho_D,r1,mL1,cL1,drx,dry); %gives lambda values of two intersection points on the line and circle centered at rc (radius 2*rho_D)
                        lambdac=lambdac0(0<lambdac0 & lambdac0<1);
                        if length(lambdac)==2
                            %lambda2=max(lambdac);
                            lam11=min(lambdac);
                        elseif length(lambdac)==1
                            lam11=lambdac;
                        elseif lambdac0(1)<0 && lambdac0(2)<0
                            angc=angInterSecCircCirc(r1,2*rho_D,rVC2,rho_safe,rotMat);  %angles of intersection points between two circles with respect to first end point
                            angfi1=max(0,min(angc));
                            lam11=0;
                        elseif max(lambdac0)>1 && min(lambdac0)<0
                            lam12=1;
                            lam11=0;
                        end
                    else
                        angfi1=angi1;
                        lam11=lambda1;
                        if lambda1<0
                            angc=angInterSecCircCirc(r1,2*rho_D,rVC2,rho_safe,rotMat);  %angles of intersection points between two circles with respect to first end point
                            angfi1=max(angi1,min(angc));
                            lam11=0;
                        end
                    end
                else
                    angc=angInterSecCircCirc(r2,2*rho_D,rVC2,rho_safe,rotMat);  %angles of intersection points between two circles with respect to first end point
                    lambdac=lambdaInterSecCircLine(r22,2*rho_D,r1,mL1,cL1,drx,dry); %gives lambda values of two intersection points on the line and circle centered at rc (radius 2*rho_D)
                    lam11=max(0,min(lambdac));
                    angfi1=max(0,min(angc));
                    angfi2=ang2;
                    lam12=1;
                end
            end
        end
    elseif angi0>0 && angi0<ang2
        if lambda0>0 && lambda0<1
            if lambda1<0
                %Check which side of the circular arc the point r1 is on
                r1_temp=rotMat*(r1-rVC2);
                angr1=atan2(r1_temp(2),r1_temp(1));
                if angr1<angi0  %on the side towards r21
                    angc=angInterSecCircCirc(r1,2*rho_D,rVC2,rho_safe,rotMat);
                    if ~isempty(angc)
                    angi1=max(angi1,min(angc));
                    end
                else
                    angc=angInterSecCircCirc(r1,2*rho_D,rVC2,rho_safe,rotMat);
                    if ~isempty(angc)
                    angi1=min(angi1,max(angc));
                    end
                end
            end
            if lambda2>1
                %Check which side of the circular arc the point r1 is on
                r2_temp=rotMat*(r2-rVC2);
                angr2=atan2(r2_temp(2),r2_temp(1));
                if angr2>angi0  %on the side towards r22
                    angc=angInterSecCircCirc(r2,2*rho_D,rVC2,rho_safe,rotMat);
                    if ~isempty(angc)
                       angi2=min(angi2,max(angc));
                    end
                else
                    angc=angInterSecCircCirc(r1,2*rho_D,rVC2,rho_safe,rotMat);
                    if ~isempty(angc)
                    angi2=max(angi2,min(angc));
                    end
                end
            end
            if angi1<0
                lambdac=lambdaInterSecCircLine(r21,2*rho_D,r1,mL1,cL1,drx,dry);
                %lambda1=max(lambda1,min(lambdac));
                if ~isempty(lambdac)
                lambda1=min(lambdac);
                end
            elseif angi1>ang2
                lambdac=lambdaInterSecCircLine(r22,2*rho_D,r1,mL1,cL1,drx,dry);
               % lambda1=max(lambda1,min(lambdac));
               if ~isempty(lambdac)
                lambda1=min(lambdac);
               end
            end
            if angi2>ang2
                lambdac=lambdaInterSecCircLine(r22,2*rho_D,r1,mL1,cL1,drx,dry);
                if ~isempty(lambdac)
                lambda2=max(lambdac);
                end
            elseif angi2<0
                lambdac=lambdaInterSecCircLine(r21,2*rho_D,r1,mL1,cL1,drx,dry);
                if ~isempty(lambdac)
                lambda2=max(lambdac);
                end
            end
            angfi1=max(0,min(angi1,angi2));
            angfi2=min(ang2,max(angi1,angi2));
            lam11=max(0,lambda1);
            lam12=min(1,lambda2);
            interSec.Flag=1;
            
        elseif lambda0<=0
            if norm(r1-rVC2)<=rho_safe+2*rho_D
                angc=angInterSecCircCirc(r1,2*rho_D,rVC2,rho_safe,rotMat);
                if (0<angc(1) && angc(1)<ang2) || (0<angc(2) && angc(2)<ang2) || norm(r1-r21)<=2*rho_D || norm(r1-r22)<=2*rho_D  %collision exists
                    interSec.Flag=1;
                    
                    r1_temp=rotMat*(r1-rVC2);
                    angr1=atan2(r1_temp(2),r1_temp(1));
                    if angr1>angi0  %Check which side of the circular arc the point r1 is on
                        lam11=0;
                        angc=angc(0<angc & angc<ang2);
                        if ~isempty(angc)
                            angfi1=max(0,min(angc));
                        else
                            angfi1=0;
                        end
                        if lambda2>1
                            angc=angInterSecCircCirc(r2,2*rho_D,rVC2,rho_safe,rotMat);
                            angi2=min(angi2,max(angc));
                        end
                        if angi2>ang2
                            lambdac0=lambdaInterSecCircLine(r22,2*rho_D,r1,mL1,cL1,drx,dry);
                            %lambda2=min(lambda2,max(lambdac));
                            lambdac=lambdac0(0<lambdac0 & lambdac0<1);
                            if ~isempty(lambdac)
                                lambda2=max(lambdac);
                            elseif max(lambdac0)>1 && min(lambdac0)<0
                                lambda2=1;
                            end
                        end
                        angfi2=min(ang2,angi2);
                        lam12=min(1,lambda2);
                    else
                        lam11=0;
                        angc=angc(0<angc & angc<ang2);
                        if ~isempty(angc)
                            angfi2=min(ang2,max(angc));
                        else
                            angfi2=ang2;
                        end
                        if lambda2>1
                            angc=angInterSecCircCirc(r2,2*rho_D,rVC2,rho_safe,rotMat);
                            angi1=max(angi1,min(angc));
                        end
                        if angi1>ang2
                            lambdac0=lambdaInterSecCircLine(r22,2*rho_D,r1,mL1,cL1,drx,dry);
                            %lambda2=min(lambda2,max(lambdac));
                            lambdac=lambdac0(0<lambdac0 & lambdac0<1);
                            if ~isempty(lambdac)
                                lambda2=max(lambdac);
                            elseif max(lambdac0)>1 && min(lambdac0)<0
                                lambda2=1;
                            end
                        end
                        angfi1=max(0,angi2);
                        lam12=min(1,lambda2);
                    end
                end
            end
        elseif lambda0>=1
            if norm(r2-rVC2)<=rho_safe+2*rho_D
                angc=angInterSecCircCirc(r2,2*rho_D,rVC2,rho_safe,rotMat);
                if (0<angc(1) && angc(1)<ang2) || (0<angc(2) && angc(2)<ang2) || norm(r2-r21)<=2*rho_D || norm(r2-r22)<=2*rho_D  %collision exists %collision exists
                    interSec.Flag=1;
                    
                    r2_temp=rotMat*(r2-rVC2);
                    angr2=atan2(r2_temp(2),r2_temp(1));
                    if angr2>angi0 %Check which side of the circular arc the point r2 is on
                        lam12=1;
                        angc=angc(0<=angc & angc<=ang2);
                        if ~isempty(angc)
                            angfi1=max(0,min(angc));
                        else
                            angfi1=0;
                        end
                        
                        if lambda1<0
                            angc=angInterSecCircCirc(r1,2*rho_D,rVC2,rho_safe,rotMat);
                            angi2=min(angi2,max(angc));
                        end
                        if angi1>ang2
                            lambdac0=lambdaInterSecCircLine(r22,2*rho_D,r1,mL1,cL1,drx,dry);
                            %lambda2=min(lambda2,max(lambdac));
                            lambdac=lambdac0(0<=lambdac0 & lambdac0<=1);
                            if ~isempty(lambdac)
                                lambda1=min(lambdac);
                            elseif max(lambdac0)>1 && min(lambdac0)<0
                                lambda1=0;
                            end
                        end
                        angfi2=min(ang2,angi1);
                        lam11=max(0,lambda1);
                    else
                        lam12=1;
                        angc=angc(0<=angc & angc<=ang2);
                        if ~isempty(angc)
                            angfi2=min(ang2,max(angc));
                        else
                            angfi2=ang2;
                        end
                        
                        if lambda1<0
                            angc=angInterSecCircCirc(r1,2*rho_D,rVC2,rho_safe,rotMat);
                            angi1=max(0,min(angc));
                        end
                        if angi2>ang2
                            lambdac0=lambdaInterSecCircLine(r22,2*rho_D,r1,mL1,cL1,drx,dry);
                            %lambda2=min(lambda2,max(lambdac));
                            lambdac=lambdac0(0<=lambdac0 & lambdac0<=1);
                            if ~isempty(lambdac)
                                lambda1=min(lambdac);
                            elseif max(lambdac0)>1 && min(lambdac0)<0
                                lambda1=0;
                            end
                        end
                        angfi1=max(0,angi1);
                        lam11=max(0,lambda1);
                    end
                end
            end
        end
    else %angi0<0
        if lambda0>0 && lambda0<1            
            rp=r21;
            if drx~=0
                x_p=(rp(1)+mL1*(rp(2)-cL1))/(1+mL1^2);
                lambdap=(x_p-r1(1))/drx;
                drpL1=((cL1+mL1*rp(1)-rp(2))/sqrt(1+mL1^2));
            else
                y_p=rp(2);
                lambdap=(y_p-r1(2))/dry;
                drpL1=abs(r1(1)-rp(1));
            end
            dlambdai=2*rho_D/L1;
            %check1=-dthetai<=angr1 && angr1<=ang2+dthetai;
            check2=(-dlambdai<=lambdap && lambdap<=1+dlambdai);
            if abs(cL1+mL1*r21(1)-r21(2))/sqrt(1+mL1^2) <2*rho_D && check2 %distance between the 2nd end point on the arc and the line
                interSec.Flag=1;
                angfi1=0;
                %Find two intersection points of the circle at the
                %end point of arc and segment
                rc=r21;
                lambdac=lambdaInterSecCircLine(rc,2*rho_D,r1,mL1,cL1,drx,dry); %gives lambda values of two intersection points on the line and circle centered at rc (radius 2*rho_D)
                
                %Check which side of the circular arc the point r1 is on
                r1_temp=rotMat*(r1-rVC2);
                angr1=atan2(r1_temp(2),r1_temp(1));
                if angr1<angi0   %
                    if lambda2>=1
                        lam12=1;
                        lam11=max(0,min(lambdac));
                        %Find two intersection points of the circle at the
                        %end point of segment 2 and segment 1
                        rc=r2;
                        angc=angInterSecCircCirc(rc,2*rho_D,rVC2,rho_safe,rotMat);  %angles of intersection points between two circles with respect to first end point
                        angfi2=min(ang2,max(angc));
                    else
                        lam11=max(0,min(lambdac));
                        if ang2<=angi2
                            angfi2=ang2;
                            lambdac=lambdaInterSecCircLine(r22,2*rho_D,r1,mL1,cL1,drx,dry); %gives lambda values of two intersection points on the line and circle centered at rc (radius 2*rho_D)
                            lam12=min(1,max(lambdac));
                        else
                            lam12=lambda2;
                            angfi2=angi2;
                        end
                    end
                else
                    if lambda1<=0
                        lam11=0;
                        lam12=min(1,max(lambdac));
                        %Find two intersection points of the circle at the
                        %end point of segment 2 and segment 1
                        rc=r1;
                        angc=angInterSecCircCirc(rc,2*rho_D,rVC2,rho_safe,rotMat);  %angles of intersection points between two circles with respect to first end point
                        angfi2=min(ang2,max(angc));
                    else
                        lam12=min(1,max(lambdac));
                        if ang2<=angi1
                            angfi2=ang2;
                            lambdac=lambdaInterSecCircLine(r22,2*rho_D,r1,mL1,cL1,drx,dry); %gives lambda values of two intersection points on the line and circle centered at rc (radius 2*rho_D)
                            lam11=max(0,min(lambdac));
                        else
                            lam11=lambda1;
                            angfi2=angi1;
                        end
                    end
                end
            end
            
        elseif lambda0<=0
            r1_temp=rotMat*(r1-rVC2);
            angr1=atan2(r1_temp(2),r1_temp(1));
            rp=r21;
            if drx~=0
                x_p=(rp(1)+mL1*(rp(2)-cL1))/(1+mL1^2);
                lambdap=(x_p-r1(1))/drx;
                drpL1=((cL1+mL1*rp(1)-rp(2))/sqrt(1+mL1^2));
            else
                y_p=rp(2);
                lambdap=(y_p-r1(2))/dry;
                drpL1=abs(r1(1)-rp(1));
            end
            dlambdai=2*rho_D/L1;
            check1=-dthetai<=angr1 && angr1<=ang2+dthetai;
            check2=(-dlambdai<=lambdap && lambdap<=1+dlambdai);
            if norm(r1-rVC2)<=rho_safe+2*rho_D && drpL1<=2*rho_D && (check1 || check2) %norm(r21-r1)<=2*rho_D   %(0<=lambdac(1)&& lambdac(1)<=1) || (0<=lambdac(2)&& lambdac(2)<=1)
                interSec.Flag=1;
                if angr1>angi0   %the segment is toward the arc
                    angfi1=0;
                    lambdac0=lambdaInterSecCircLine(r21,2*rho_D,r1,mL1,cL1,drx,dry); %gives lambda values of two intersection points on the line and circle centered at rc (radius 2*rho_D)
                    lambdac=lambdac0(0<lambdac0 & lambdac0<1);
                    if length(lambdac)==2
                        %lambda2=max(lambdac);
                        lam11=min(lambdac);
                    elseif length(lambdac)==1
                        lam11=0;
                    elseif lambdac0(1)<0 && lambdac0(2)<0
                        angc=angInterSecCircCirc(r1,2*rho_D,rVC2,rho_safe,rotMat);  %angles of intersection points between two circles with respect to first end point
                        angfi1=max(0,min(angc));
                        lam11=0;
                    elseif max(lambdac0)>1 && min(lambdac0)<0
                        lam12=1;
                        lam11=0;
                    end
                    
                    if angi2>ang2
                        angfi2=ang2;
                        lambdac0=lambdaInterSecCircLine(r22,2*rho_D,r1,mL1,cL1,drx,dry); %gives lambda values of two intersection points on the line and circle centered at rc (radius 2*rho_D)
                        lambdac=lambdac0(0<lambdac0 & lambdac0<1);
                        if length(lambdac)==2
                            %lambda2=max(lambdac);
                            lam12=max(lambdac);
                        elseif length(lambdac)==1
                            lam12=lambdac;
                        elseif lambdac0(1)>1 && lambdac0(2)>1
                            angc=angInterSecCircCirc(r2,2*rho_D,rVC2,rho_safe,rotMat);  %angles of intersection points between two circles with respect to first end point
                            angfi2=min(ang2,max(angc));
                            lam12=1;
                        elseif max(lambdac0)>1 && min(lambdac0)<0
                            lam12=1;
                            lam11=0;
                        end
                    else
                        angfi2=angi2;
                        lam12=lambda2;
                        if lambda2>1
                            angc=angInterSecCircCirc(r2,2*rho_D,rVC2,rho_safe,rotMat);  %angles of intersection points between two circles with respect to first end point
                            angfi2=min(angi2,max(angc));
                            lam12=1;
                        end
                    end
                else
                    angc=angInterSecCircCirc(r1,2*rho_D,rVC2,rho_safe,rotMat);  %angles of intersection points between two circles with respect to first end point
                    lambdac=lambdaInterSecCircLine(r21,2*rho_D,r1,mL1,cL1,drx,dry); %gives lambda values of two intersection points on the line and circle centered at rc (radius 2*rho_D)
                    lam12=min(1,max(lambdac));
                    angfi2=min(ang2,max(angc));
                    angfi1=0;
                    lam11=0;
                end
            end
        elseif lambda0>=1
            r2_temp=rotMat*(r2-rVC2);
            angr2=atan2(r2_temp(2),r2_temp(1));
            rp=r21;
            if drx~=0
                x_p=(rp(1)+mL1*(rp(2)-cL1))/(1+mL1^2);
                lambdap=(x_p-r1(1))/drx;
                drpL1=((cL1+mL1*rp(1)-rp(2))/sqrt(1+mL1^2));
            else
                y_p=rp(2);
                lambdap=(y_p-r1(2))/dry;
                drpL1=abs(r1(1)-rp(1));
            end
            dlambdai=2*rho_D/L1;
            check1=-dthetai<=angr2 && angr2<=ang2+dthetai;
            check2=(-dlambdai<=lambdap && lambdap<=1+dlambdai);
            if norm(r2-rVC2)<=rho_safe+2*rho_D && drpL1<=2*rho_D && (check1 || check2) %norm(r21-r2)<=2*rho_D
                interSec.Flag=1;
                if angr2>angi0
                    angfi1=0;
                    lambdac0=lambdaInterSecCircLine(r21,2*rho_D,r1,mL1,cL1,drx,dry); %gives lambda values of two intersection points on the line and circle centered at rc (radius 2*rho_D)
                    lambdac=lambdac0(0<lambdac0 & lambdac0<1);
                    if length(lambdac)==2
                        %lambda2=max(lambdac);
                        lam12=max(lambdac);
                    elseif length(lambdac)==1
                        lam12=1;
                    elseif lambdac0(1)>1 && lambdac0(2)>1
                        angc=angInterSecCircCirc(r2,2*rho_D,rVC2,rho_safe,rotMat);  %angles of intersection points between two circles with respect to first end point
                        angfi1=max(0,min(angc));
                        lam12=1;
                    elseif max(lambdac0)>1 && min(lambdac0)<0
                        lam12=1;
                        lam11=0;
                    end
                    
                    if angi1>ang2
                        angfi2=ang2;
                        lambdac0=lambdaInterSecCircLine(r22,2*rho_D,r1,mL1,cL1,drx,dry); %gives lambda values of two intersection points on the line and circle centered at rc (radius 2*rho_D)
                        lambdac=lambdac0(0<lambdac0 & lambdac0<1);
                        if length(lambdac)==2
                            %lambda2=max(lambdac);
                            lam11=min(lambdac);
                        elseif length(lambdac)==1
                            lam11=lambdac;
                        elseif lambdac0(1)<0 && lambdac0(2)<0
                            angc=angInterSecCircCirc(r1,2*rho_D,rVC2,rho_safe,rotMat);  %angles of intersection points between two circles with respect to first end point
                            angfi2=min(ang2,max(angc));
                            lam11=0;
                        elseif max(lambdac0)>1 && min(lambdac0)<0
                            lam12=1;
                            lam11=0;
                        end
                    else
                        angfi2=angi1;
                        lam11=lambda1;
                        if lambda1<0
                            angc=angInterSecCircCirc(r1,2*rho_D,rVC2,rho_safe,rotMat);  %angles of intersection points between two circles with respect to first end point
                            angfi2=min(angi1,max(angc));
                            lam11=0;
                        end
                    end
                else
                    angc=angInterSecCircCirc(r2,2*rho_D,rVC2,rho_safe,rotMat);  %angles of intersection points between two circles with respect to first end point
                    lambdac=lambdaInterSecCircLine(r21,2*rho_D,r1,mL1,cL1,drx,dry); %gives lambda values of two intersection points on the line and circle centered at rc (radius 2*rho_D)
                    lam11=max(0,min(lambdac));
                    angfi2=min(ang2,max(angc));
                    angfi1=0;
                    lam12=1;
                end
            end
        end
    end
    
    
    %                 interSec.Pos2(:,1)=rVC2(:)+(rho_safe)*[cos(angfi1+ang1),sin(angfi1+ang1)]'; %first intersection point on path 2
    %                 interSec.Pos2(:,2)=rVC2(:)+(rho_safe)*[cos(angfi2+ang1),sin(angfi2+ang1)]'; %second intersection point on path 2
    if interSec.Flag==1
        if ~exist('lam11','var')
            1;
        end
        
        interSec.Pos1(:,1)=[r1(1)+lam11*drx, r1(2)+lam11*dry]'; %first intersection point on path 1
        interSec.Pos1(:,2)=[r1(1)+lam12*drx, r1(2)+lam12*dry]'; %second intersection point on path 1
        
        ds1=norm(interSec.Pos1(:,1)-interSec.Pos1(:,2));
        interSec.Pbar1=ds1;
        interSec.S10(1)=norm(interSec.Pos1(:,1)-r1);
        interSec.S10(2)=interSec.S10(1)+ds1;       
        
        if segType2==1            
            interSec.Pos2(:,1)=rVC2(:)+(rho_safe)*[cos(angfi1+ang1),sin(angfi1+ang1)]'; %first intersection point on path 2
            interSec.Pos2(:,2)=rVC2(:)+(rho_safe)*[cos(angfi2+ang1),sin(angfi2+ang1)]'; %second intersection point on path 2
            ds2=rho_safe*acos(1-norm(interSec.Pos2(:,1)-interSec.Pos2(:,2))^2/(2*rho_safe^2));
            interSec.Pbar2=ds2;
            interSec.S20(1)=rho_safe*acos(1-norm(interSec.Pos2(:,1)-r21)^2/(2*rho_safe^2));
            interSec.S20(2)=interSec.S20(1)+ds2;
        else
            interSec.Pos2(:,1)=rVC2(:)+(rho_safe)*[cos(angfi2+ang1),sin(angfi2+ang1)]'; %first intersection point on path 2
            interSec.Pos2(:,2)=rVC2(:)+(rho_safe)*[cos(angfi1+ang1),sin(angfi1+ang1)]'; %second intersection point on path 2
            ds2=rho_safe*acos(1-norm(interSec.Pos2(:,1)-interSec.Pos2(:,2))^2/(2*rho_safe^2));
            interSec.Pbar2=ds2;
            interSec.S20(1)=rho_safe*acos(1-norm(interSec.Pos2(:,1)-r22)^2/(2*rho_safe^2));   %r22 is the original first point
            interSec.S20(2)=interSec.S20(1)+ds2;
        end
        return;
    end
    %    end
end

end