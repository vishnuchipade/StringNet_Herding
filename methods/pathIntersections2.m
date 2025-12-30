function [interSec]=pathIntersections2(Path1,Path2)
% ------------------------------
% Author: Vishnu S. Chipade
% Email: vishnuc@umich.edu
% ------------------------------

%This function returns the intersections of the two paths when a robot of
%finite size (rho_D) moves on them

global rho_safe rho_D

NS1=length(Path1.P);
NS2=length(Path2.P);
interSec.Flag=zeros(NS1,NS2);
interSec.Pbar1=0;  %Length of the path 1 in conflict
interSec.Pbar2=0;  %Length of the path 2 in conflict
rV1=Path1.rV;
rV2=Path2.rV;
rVC1=Path1.rVC;
rVC2=Path2.rVC;
segType1=Path1.segType;
segType2=Path2.segType;

interSec.flag0=0;
%check for intersection between straight segments of first path with the
%segments of the second
for i=1:2:NS1
    
    Si0=Path1.S(i);
    
    r1=rV1(1:2,i);
    r2=rV1(1:2,i+1);
    dr=r2-r1;
    L1=norm(dr);   %length of the segment
    drx=dr(1);
    dry=dr(2);
    theta1=atan2(dry,drx);
    if theta1<0
        theta1=theta1+2*pi;
    end
    mL1=tan(theta1);
    cL1=r1(2)-mL1*r1(1);
    dr_hat=dr/norm(dr);
    
    %with straight segments of the second path
    for j=1:2:NS2
        
        Sj0=Path2.S(j);
        
        r21=rV2(1:2,j);
        r22=rV2(1:2,j+1);
        drx2=(r22(1)-r21(1));
        dry2=r22(2)-r21(2);
        theta2=atan2(dry2, drx2);
        if theta2<0
            theta2=theta2+2*pi;
        end
        dtheta=theta2-theta1;
        if dtheta<0
            dtheta=dtheta+2*pi;
        end
        mL2=tan(theta2);
        cL2=r21(2)-mL2*r21(1);
        L2=norm(r21-r22); %length of the segment
        if mL1==mL2
            %&& (r1(1)==r21(1) && r1(2)==r21(2)) || (r1(1)==r22(1) && r1(2)==r22(2))
            %find the distance between the two || lines
            d=abs(cL1-cL2)/sqrt(1+mL1^2);
            if d<2*rho_D  %intersection of the rectangular area swept by the footprint of the defender
                interSec.Flag(i,j)=2;   %2 if the intersection is happening over a section or entire segment
                if L1<L2
                    interSec.Pos1{i,j}(:,1)=r1;  %first end point of the intersection segment on path 1
                    interSec.Pos1{i,j}(:,2)=r2;   %second end point of the intersection segment on path 1
                    interSec.S10{i,j}(1)=Si0+0;
                    interSec.S10{i,j}(2)=Si0+L1;
                    interSec.Pbar1=interSec.Pbar1+L1;
                    
                    interSec.Pos2{i,j}(1,1)=(r1(1)+mL1*r1(2)-mL1*cL2)/(1+mL1^2);  %x-coord of first end point of the intersection segment on path 2
                    interSec.Pos2{i,j}(1,2)=(r2(1)+mL1*r2(2)-mL1*cL2)/(1+mL1^2);  %x-coord of second end point of the intersection segment on path 2
                    if mL2<1e16
                        interSec.Pos2{i,j}(2,1)=mL2*interSec.Pos2{i,j}(1,1)+cL2;
                        interSec.Pos2{i,j}(2,2)=mL2*interSec.Pos2{i,j}(1,2)+cL2;
                    else
                        interSec.Pos2{i,j}(2,1)=r1(2);
                        interSec.Pos2{i,j}(2,2)=r2(2);
                    end
                    interSec.S20{i,j}(1)=Sj0+norm(r21-interSec.Pos2{i,j}(:,1));
                    interSec.S20{i,j}(2)=interSec.S20{i,j}(1)+L1;
                    interSec.Pbar2=interSec.Pbar2+L1;
                else
                    interSec.Pos2{i,j}(:,1)=r21;  %first end point of the intersection segment on path 2
                    interSec.Pos2{i,j}(:,2)=r22;   %second end point of the intersection segment on path 2
                    interSec.S20{i,j}(1)=Sj0+0;
                    interSec.S20{i,j}(2)=Sj0+L2;
                    interSec.Pbar2=interSec.Pbar2+L2;
                    
                    interSec.Pos1{i,j}(1,1)=(r21(1)+mL1*r21(2)-mL1*cL2)/(1+mL1^2);  %x-coord of  first end point of the intersection segment on path 1
                    interSec.Pos1{i,j}(1,2)=(r22(1)+mL1*r22(2)-mL1*cL2)/(1+mL1^2);  %x-coord of  second end point of the intersection segment on path 1
                    if mL1<1e16
                        interSec.Pos1{i,j}(2,1)=mL1*interSec.Pos1{i,j}(1,1)+cL1;
                        interSec.Pos1{i,j}(2,2)=mL1*interSec.Pos1{i,j}(1,2)+cL1;
                    else
                        interSec.Pos1{i,j}(2,1)=r21(2);
                        interSec.Pos1{i,j}(2,2)=r22(2);
                    end
                    interSec.S10{i,j}(1)=Si0+norm(r1-interSec.Pos1{i,j}(:,1));
                    interSec.S10{i,j}(2)=interSec.S10{i,j}(1)+L2;
                    interSec.Pbar1=interSec.Pbar1+L2;
                end
                break;
            end
        else
            %cL2=r21(2)-mL2*r21(1);
            x_int=(cL2-cL1)/(mL1-mL2);
            if mL1>1e16
                y_int=mL2*x_int+cL2;
            else
                y_int=mL1*x_int+cL1;
            end
            Pbar0=max(rho_D*sqrt(2/(1-abs(cos(dtheta)))),2*rho_D);
            dx=Pbar0*abs(cos(theta1));
            dy=Pbar0*abs(sin(theta1));
            xbar11=x_int-sign(drx)*dx;
            ybar11=y_int-sign(dry)*dy;
            xbar12=x_int+sign(drx)*dx;
            ybar12=y_int+sign(dry)*dy;
            
            dx=Pbar0*abs(cos(theta2));
            dy=Pbar0*abs(sin(theta2));
            xbar21=x_int-sign(drx2)*dx;
            ybar21=y_int-sign(dry2)*dy;
            xbar22=x_int+sign(drx2)*dx;
            ybar22=y_int+sign(dry2)*dy;
            
            %y_int=mL1*x_int+cL1;
            if drx~=0
                lambda1=(x_int-r1(1))/drx;%norm([x_int,y_int]'-r1)/L1;
                lambda11=(xbar11-r1(1))/drx;%norm([xbar11,ybar11]'-r1)/L1;
                lambda12=(xbar12-r1(1))/drx;%norm([xbar12,ybar12]'-r1)/L1;
            else
                lambda1=(y_int-r1(2))/dry;%norm([x_int,y_int]'-r1)/L1;
                lambda11=(ybar11-r1(2))/dry;%norm([xbar11,ybar11]'-r1)/L1;
                lambda12=(ybar12-r1(2))/dry;%norm([xbar12,ybar12]'-r1)/L1;
            end
            if drx2~=0
                lambda2=(x_int-r21(1))/drx2;%norm([x_int,y_int]'-r21)/L2;
                lambda21=(xbar21-r21(1))/drx2;%norm([xbar21,ybar21]'-r21)/L2;
                lambda22=(xbar22-r21(1))/drx2;%norm([xbar22,ybar22]'-r21)/L2;
            else
                lambda2=(y_int-r21(2))/dry2;%norm([x_int,y_int]'-r21)/L2;
                lambda21=(ybar21-r21(2))/dry2;%norm([xbar21,ybar21]'-r21)/L2;
                lambda22=(ybar22-r21(2))/dry2;%norm([xbar22,ybar22]'-r21)/L2;
            end
            %flagIntersect=((0<= lambda11) && (1>= lambda11)) || ((0<= lambda12) && (1>= lambda12)) || ...
            %((0<= lambda21) && (1>= lambda21)) ||  ((0<= lambda22) && (1>= lambda22));
            
            flagIntersect=0;
            %             if (0<= lambda1) && (0<= lambda2)   %both path's end part are possibly in conflict
            %                 if (lambda11<=1 || lambda12<=1) && (lambda21<=1 || lambda22<=1)
            %                     flagIntersect=1;
            %                 end
            %             elseif (0<= lambda1) && (0> lambda2)
            %                 if (lambda11<=1 || lambda12<=1) && (lambda21>0 || lambda22>=0)
            %                     flagIntersect=1;
            %                 end
            %             elseif (0> lambda1) && (0<= lambda2)
            %                 if (lambda11>0 || lambda12>=0) && (lambda21<=1 || lambda22<=1)
            %                     flagIntersect=1;
            %                 end
            %             else %(0> lambda1) && (0> lambda2)
            %                 if (lambda11>0 || lambda12>=0) && (lambda21>0 || lambda22>=0)
            %                     flagIntersect=1;
            %                 end
            %             end
            if lambda2>=1
                if lambda1>0 && lambda1<1
                    % interSec.Flag(i,j)=1;   %1 if point intersection
                    if norm(r22-[x_int,y_int]')<2*rho_D/abs(sin(dtheta))
                        lam22=1;
                        %Find two intersection points of the circle at the
                        %end point of segment 2 and segment 1
                        rc=r22;
                        lambda0=lambdaInterSecCircLine(rc,2*rho_D,r1,mL1,cL1,drx,dry); %gives lambda values of two intersection points on the line and circle centered at rc (radius 2*rho_D)
                        if (r21-[x_int,y_int]')'*(r2-[x_int,y_int]')>0
                            lam11=max(0,min(lambda0));
                            if lambda12>1
                                lam12=1;
                                %Find two intersection points of the circle at the
                                %end point of segment 2 and segment 1
                                rc=r2;
                                lambda00=lambdaInterSecCircLine(rc,2*rho_D,r21,mL2,cL2,drx2,dry2); %gives lambda values of two intersection points on the line and circle centered at rc (radius 2*rho_D)
                                lam21=max(0,min(lambda00));
                            else
                                lam12=lambda12;
                                lam21=max(0,lambda21);
                            end
                        else
                            lam12=min(1,max(lambda0));
                            if lambda11<0
                                lam11=0;
                                %Find two intersection points of the circle at the
                                %end point of segment 2 and segment 1
                                rc=r1;
                                lambda00=lambdaInterSecCircLine(rc,2*rho_D,r21,mL2,cL2,drx2,dry2); %gives lambda values of two intersection points on the line and circle centered at rc (radius 2*rho_D)
                                lam21=max(0,min(lambda00));
                            else
                                lam11=lambda11;
                                lam21=max(0,lambda21);
                            end
                        end
                        flagIntersect=1;
                    end
                elseif lambda1<=0
                    if norm(r22-r1)<2*rho_D
                        lam11=0;
                        lam22=1;
                        if (r21-[x_int,y_int]')'*(r2-[x_int,y_int]')>0
                            lam12=min(1,lambda12);
                            lam21=max(0,lambda21);
                        else
                            a=norm(r22-[x_int,y_int]');
                            %b=a*abs(cos(dtheta))+sqrt(4*(rho_D)^2-a^2*sin(dtheta)^2);
                            b=-a*cos(dtheta)+sqrt(4*(rho_D)^2-a^2*sin(dtheta)^2);
                            lam12=min(1,lambda1+b/L1);
                            a=norm(r1-[x_int,y_int]');
                            %b=a*abs(cos(dtheta))+sqrt(4*(rho_D)^2-a^2*sin(dtheta)^2);
                            b=-a*cos(dtheta)+sqrt(4*(rho_D)^2-a^2*sin(dtheta)^2);
                            lam21=max(0,lambda2-b/L2);
                        end
                        
                        flagIntersect=1;
                    end
                elseif lambda1>=1
                    if norm(r22-r2)<2*rho_D
                        lam12=1;
                        lam22=1;
                        if (r21-[x_int,y_int]')'*(r1-[x_int,y_int]')>0
                            lam11=max(lambda11,0);
                            lam21=max(lambda21,0);
                        else
                            a=norm(r22-[x_int,y_int]');
                            %b=a*abs(cos(dtheta))+sqrt(4*(rho_D)^2-a^2*sin(dtheta)^2);
                            b=-a*cos(dtheta)+sqrt(4*(rho_D)^2-a^2*sin(dtheta)^2);
                            lam11=max(1,lambda1-b/L1);
                            a=norm(r2-[x_int,y_int]');
                            %b=a*abs(cos(dtheta))+sqrt(4*(rho_D)^2-a^2*sin(dtheta)^2);
                            b=-a*cos(dtheta)+sqrt(4*(rho_D)^2-a^2*sin(dtheta)^2);
                            lam21=max(0,lambda2-b/L2);
                        end
                        
                        flagIntersect=1;
                    end
                end
            elseif lambda2>0 && lambda2<1
                if lambda1>0 && lambda1<1
                    %interSec.Flag(i,j)=1;   %1 if point intersection
                    if lambda11<0
                        a=norm(r1-[x_int,y_int]');
                        %b=a*abs(cos(dtheta))+sqrt(4*(rho_D)^2-a^2*sin(dtheta)^2);
                        b=a*cos(dtheta)+sqrt(4*(rho_D)^2-a^2*sin(dtheta)^2);  %positve cos(dtheta) term
                        if (r1-[x_int,y_int]')'*(r21-[x_int,y_int]')>0
                            lambda21=max(lambda21,lambda2-b/L2);
                        else
                            lambda22=min(lambda22,lambda2+b/L2);
                        end
                    end
                    if lambda12>1
                        a=norm(r2-[x_int,y_int]');
                        %b=a*abs(cos(dtheta))+sqrt(4*(rho_D)^2-a^2*sin(dtheta)^2);
                        b=a*cos(dtheta)+sqrt(4*(rho_D)^2-a^2*sin(dtheta)^2);
                        if (r2-[x_int,y_int]')'*(r21-[x_int,y_int]')>0
                            lambda21=max(lambda21,lambda2-b/L2);
                        else
                            lambda22=min(lambda22,lambda2+b/L2);
                        end
                    end
                    if lambda21<0
                        a=norm(r21-[x_int,y_int]');
                        %b=a*abs(cos(dtheta))+sqrt(4*(rho_D)^2-a^2*sin(dtheta)^2);
                        b=a*cos(dtheta)+sqrt(4*(rho_D)^2-a^2*sin(dtheta)^2);
                        if (r21-[x_int,y_int]')'*(r1-[x_int,y_int]')>0
                            lambda11=max(lambda11,lambda1-b/L1);
                        else
                            lambda12=min(lambda12,lambda1+b/L1);
                        end
                    end
                    if lambda22>1
                        a=norm(r22-[x_int,y_int]');
                        %b=a*abs(cos(dtheta))+sqrt(4*(rho_D)^2-a^2*sin(dtheta)^2);
                        b=a*cos(dtheta)+sqrt(4*(rho_D)^2-a^2*sin(dtheta)^2);
                        if (r22-[x_int,y_int]')'*(r1-[x_int,y_int]')>0
                            lambda11=max(lambda11,lambda1-b/L1);
                        else
                            lambda12=min(lambda12,lambda1+b/L1);
                        end
                    end
                    lam11=max(lambda11,0);
                    lam12=min(1,lambda12);
                    lam21=max(lambda21,0);
                    lam22=min(1,lambda22);
                    flagIntersect=1;
                elseif lambda1<=0
                    if norm(r1-[x_int,y_int]')<2*rho_D/abs(sin(dtheta))
                        lam11=0;
                        %Find two intersection points of the circle at the
                        %end point of segment 2 and segment 1
                        rc=r1;
                        lambda0=lambdaInterSecCircLine(rc,2*rho_D,r21,mL2,cL2,drx2,dry2); %gives lambda values of two intersection points on the line and circle centered at rc (radius 2*rho_D)
                        if (r22-[x_int,y_int]')'*(r2-[x_int,y_int]')>0
                            lam21=max(0,min(lambda0));
                            if lambda22>1
                                lam22=1;
                                %Find two intersection points of the circle at the
                                %end point of segment 2 and segment 1
                                rc=r22;
                                lambda00=lambdaInterSecCircLine(rc,2*rho_D,r1,mL1,cL1,drx,dry); %gives lambda values of two intersection points on the line and circle centered at rc (radius 2*rho_D)
                                lam12=min(1,max(lambda00));
                            else
                                lam22=lambda22;
                                lam12=min(1,lambda12);
                            end
                        else
                            lam22=min(1,max(lambda0));
                            if lambda21<0
                                lam21=0;
                                %Find two intersection points of the circle at the
                                %end point of segment 2 and segment 1
                                rc=r21;
                                lambda00=lambdaInterSecCircLine(rc,2*rho_D,r1,mL1,cL1,drx,dry); %gives lambda values of two intersection points on the line and circle centered at rc (radius 2*rho_D)
                                lam12=min(1,max(lambda00));
                            else
                                lam21=lambda21;
                                lam12=min(1,lambda12);
                            end
                        end
                        flagIntersect=1;
                    end
                elseif lambda1>=1
                    if norm(r2-[x_int,y_int]')<2*rho_D/abs(sin(dtheta))
                        lam12=1;
                        %Find two intersection points of the circle at the
                        %end point of segment 2 and segment 1
                        rc=r2;
                        lambda0=lambdaInterSecCircLine(rc,2*rho_D,r21,mL2,cL2,drx2,dry2); %gives lambda values of two intersection points on the line and circle centered at rc (radius 2*rho_D)
                        if (r22-[x_int,y_int]')'*(r1-[x_int,y_int]')>0
                            lam21=max(0,min(lambda0));
                            if lambda22>1
                                lam22=1;
                                %Find two intersection points of the circle at the
                                %end point of segment 2 and segment 1
                                rc=r22;
                                lambda00=lambdaInterSecCircLine(rc,2*rho_D,r1,mL1,cL1,drx,dry); %gives lambda values of two intersection points on the line and circle centered at rc (radius 2*rho_D)
                                lam11=max(0,min(lambda00));
                            else
                                lam22=lambda22;
                                lam11=max(0,lambda11);
                            end
                        else
                            lam22=min(1,max(lambda0));
                            if lambda21<0
                                lam21=0;
                                %Find two intersection points of the circle at the
                                %end point of segment 2 and segment 1
                                rc=r21;
                                lambda00=lambdaInterSecCircLine(rc,2*rho_D,r1,mL1,cL1,drx,dry); %gives lambda values of two intersection points on the line and circle centered at rc (radius 2*rho_D)
                                lam11=max(0,min(lambda00));
                            else
                                lam21=lambda21;
                                lam11=max(0,lambda11);
                            end
                        end
                        flagIntersect=1;
                    end
                end
            else %lambda2<0
                if lambda1>0 && lambda1<1
                    %interSec.Flag(i,j)=1;   %1 if point intersection
                    if norm(r21-[x_int,y_int]')<2*rho_D/abs(sin(dtheta))
                        lam21=0;
                        %Find two intersection points of the circle at the
                        %end point of segment 2 and segment 1
                        rc=r21;
                        lambda0=lambdaInterSecCircLine(rc,2*rho_D,r1,mL1,cL1,drx,dry); %gives lambda values of two intersection points on the line and circle centered at rc (radius 2*rho_D)
                        if (r22-[x_int,y_int]')'*(r2-[x_int,y_int]')>0
                            lam11=max(0,min(lambda0));
                            if lambda12>1
                                lam12=1;
                                %Find two intersection points of the circle at the
                                %end point of segment 2 and segment 1
                                rc=r2;
                                lambda00=lambdaInterSecCircLine(rc,2*rho_D,r21,mL2,cL2,drx2,dry2); %gives lambda values of two intersection points on the line and circle centered at rc (radius 2*rho_D)
                                lam22=min(1,max(lambda00));
                            else
                                lam12=lambda12;
                                lam22=min(1,lambda22);
                            end
                        else
                            lam12=min(1,max(lambda0));
                            if lambda11<0
                                lam11=0;
                                %Find two intersection points of the circle at the
                                %end point of segment 2 and segment 1
                                rc=r1;
                                lambda00=lambdaInterSecCircLine(rc,2*rho_D,r21,mL2,cL2,drx2,dry2); %gives lambda values of two intersection points on the line and circle centered at rc (radius 2*rho_D)
                                lam22=min(1,max(lambda00));
                            else
                                lam11=lambda11;
                                lam22=min(1,lambda22);
                            end
                        end
                        flagIntersect=1;
                    end
                elseif lambda1<=0
                    if norm(r21-r1)<2*rho_D
                        lam21=0;
                        lam11=0;
                        if (r22-[x_int,y_int]')'*(r1-[x_int,y_int]')>0
                            lam11=max(0,lambda11);
                            lam22=min(1,lambda22);
                        else
                            a=norm(r1-[x_int,y_int]');
                            %b=a*abs(cos(dtheta))+sqrt(4*(rho_D)^2-a^2*sin(dtheta)^2);
                            b=-a*cos(dtheta)+sqrt(4*(rho_D)^2-a^2*sin(dtheta)^2);
                            lam22=min(1,lambda2+b/L2);
                            a=norm(r21-[x_int,y_int]');
                            %b=a*abs(cos(dtheta))+sqrt(4*(rho_D)^2-a^2*sin(dtheta)^2);
                            b=-a*cos(dtheta)+sqrt(4*(rho_D)^2-a^2*sin(dtheta)^2);
                            lam12=min(1,lambda1+b/L1);
                        end
                        flagIntersect=1;
                    end
                elseif lambda1>=1
                    if norm(r2-r21)<2*rho_D
                        lam21=0;
                        lam12=1;
                        if (r22-[x_int,y_int]')'*(r1-[x_int,y_int]')>0
                            lam11=max(0,lambda11);
                            lam22=min(1,lambda22);
                        else
                            a=norm(r2-[x_int,y_int]');
                            %b=a*abs(cos(dtheta))+sqrt(4*(rho_D)^2-a^2*sin(dtheta)^2);
                            b=-a*cos(dtheta)+sqrt(4*(rho_D)^2-a^2*sin(dtheta)^2);
                            lam22=min(1,lambda2+b/L2);
                            a=norm(r21-[x_int,y_int]');
                            %b=a*abs(cos(dtheta))+sqrt(4*(rho_D)^2-a^2*sin(dtheta)^2);
                            b=-a*cos(dtheta)+sqrt(4*(rho_D)^2-a^2*sin(dtheta)^2);
                            lam11=max(0,lambda1-b/L1);
                        end
                        flagIntersect=1;
                    end
                end
            end
            
            if flagIntersect==1
                interSec.Flag(i,j)=1;
                interSec.Pos1{i,j}(:,1)=[r1(1)+lam11*drx, r1(2)+lam11*dry]';
                interSec.Pos1{i,j}(:,2)=[r1(1)+lam12*drx, r1(2)+lam12*dry]';
                interSec.Pbar1=interSec.Pbar1+(lam12-lam11)*L1;
                interSec.S10{i,j}(1)=Si0+lam11*L1;
                interSec.S10{i,j}(2)=Si0+lam12*L1;
                
                interSec.Pos2{i,j}(:,1)=[r21(1)+lam21*drx2, r21(2)+lam21*dry2]';
                interSec.Pos2{i,j}(:,2)=[r21(1)+lam22*drx2, r21(2)+lam22*dry2]';
                interSec.Pbar2=interSec.Pbar2+(lam22-lam21)*L2;
                interSec.S20{i,j}(1)=Sj0+lam21*L2;
                interSec.S20{i,j}(2)=Sj0+lam22*L2;
            end
            
        end
    end
    
    %with circular segments of the second path
    for j=2:2:NS2
        Sj0=Path2.S(j);
        r21=rV2(:,j);
        r22=rV2(:,j+1);
        interSec0=interSecLineCirc(r1,r2,mL1,cL1,dr_hat,drx,dry,r21,r22,rVC2(:,j),segType2(j));
        if interSec0.Flag~=0
            interSec.Flag(i,j)=interSec0.Flag;
            interSec.Pos1{i,j}=interSec0.Pos1;
            interSec.Pos2{i,j}=interSec0.Pos2;
            interSec.Pbar1=interSec.Pbar1+interSec0.Pbar1;
            interSec.Pbar2=interSec.Pbar1+interSec0.Pbar2;
            interSec.S10{i,j}=Si0+interSec0.S10;
            interSec.S20{i,j}=Sj0+interSec0.S20;
        end
    end
end


%check for intersection between circular segments of first path with the
%straight segments of the second
for j=1:2:NS2
    Sj0=Path2.S(j);
    r1=rV2(1:2,j);
    r2=rV2(1:2,j+1);
    dr=r2-r1;
    L1=norm(dr);   %length of the segment
    drx=r2(1)-r1(1);
    dry=r2(2)-r1(2);
    theta1=atan2(dry,drx);
    if theta1<0
        theta1=theta1+2*pi;
    end
    mL1=tan(theta1);
    cL1=r1(2)-mL1*r1(1);
    dr_hat=dr/norm(dr);
    
    for i=2:2:NS1
        Si0=Path1.S(i);
        r21=rV1(:,i);
        r22=rV1(:,i+1);
        %interSec0=interSecStraightCirc(r1,r2,mL1,cL1,dr_hat,drx,dry,r21,r22,rVC1(:,i),segType1(i));
        interSec0=interSecLineCirc(r1,r2,mL1,cL1,dr_hat,drx,dry,r21,r22,rVC1(:,i),segType1(i));
        if interSec0.Flag~=0
            interSec.Flag(i,j)=interSec0.Flag;
            interSec.Pos1{i,j}=interSec0.Pos2;  %The second one gives values for circular arc.
            interSec.Pos2{i,j}=interSec0.Pos1;
            interSec.Pbar1=interSec.Pbar1+interSec0.Pbar2;
            interSec.Pbar2=interSec.Pbar1+interSec0.Pbar1;
            interSec.S10{i,j}=Si0+interSec0.S20;
            interSec.S20{i,j}=Sj0+interSec0.S10;
        end
    end
    
end

%check for intersection between circular segments of first path with the
%circular segments of the second (Interestingly, both the arcs have to be on the same circle for them to have collision)
for i=2:2:NS1
    Si0=Path1.S(i);
    r1=rV1(:,i);
    r2=rV1(:,i+1);
    ang11=atan2(r1(2)-rVC1(2,i),r1(1)-rVC1(1,i));
    if ang11<0
        ang11=ang11+2*pi;
    end
    rotMat=[cos(ang11), sin(ang11);  -sin(ang11), cos(ang11)];
    r2_prime=rotMat*(r2-rVC1(1:2,i));   %Rotate relative to the center of arc
    ang12_prime=atan2(r2_prime(2),r2_prime(1));
    %with circular segments of the second path
    for j=2:2:NS2
        r21=rV2(:,j);
        r22=rV2(:,j+1);
        if norm(rVC1(:,i)-rVC2(:,j))<1e-5 %If the centers of the arcs are same
            %Find the final angles corresponding to the segments in
            %collision
            interSec0=interSecCircCirc(r1,r2,rotMat,ang12_prime,segType1(i),r21,r22,segType2(j),rVC1(:,i));
            angfi=interSec0.angfi;
            %ang=sort([ang11,ang12,ang21,ang22]);
            if ~isempty(angfi)
                interSec.Flag(i,j)=1;
                interSec.Pos1{i,j}(:,1)=rVC1(:,i)+rho_safe*[cos(ang11+angfi(1)), sin(ang11+angfi(1))]';  %middle angles will give the end points of the common segment
                interSec.Pos1{i,j}(:,2)=rVC1(:,i)+rho_safe*[cos(ang11+angfi(2)), sin(ang11+angfi(2))]';
                interSec.S10{i,j}(1)=Si0+rho_safe*abs(angfi(1));
                interSec.S10{i,j}(2)=Si0+rho_safe*abs(angfi(2));
                interSec.Pbar1=interSec.Pbar1+rho_safe*abs(angfi(2)-angfi(1));
                
                interSec.Pos2{i,j}(:,1)=rVC1(:,i)+rho_safe*[cos(ang11+angfi(3)), sin(ang11+angfi(3))]';  %middle angles will give the end points of the common segment
                interSec.Pos2{i,j}(:,2)=rVC1(:,i)+rho_safe*[cos(ang11+angfi(4)), sin(ang11+angfi(4))]';
                interSec.S20{i,j}(1)=Sj0+rho_safe*abs(angfi(3));
                interSec.S20{i,j}(2)=Sj0+rho_safe*abs(angfi(4));
                interSec.Pbar2=interSec.Pbar2+rho_safe*abs(angfi(4)-angfi(3));
                break;
            end
        end
    end
end

%If there are contiguous common segments on the paths then combined them
%together to treat them as a single collision segment


%Store the intersection segments together
countS1=0;
%countS2=1;
for i=1:NS1
    for j=1:NS2
        if interSec.Flag(i,j)~=0
            countS1=countS1+1;
            interSec.rV11(:,countS1)=interSec.Pos1{i,j}(:,1);
            interSec.rV12(:,countS1)=interSec.Pos1{i,j}(:,2);
            interSec.rV21(:,countS1)=interSec.Pos2{i,j}(:,1);
            interSec.rV22(:,countS1)=interSec.Pos2{i,j}(:,2);
            interSec.S11(countS1)=interSec.S10{i,j}(1);
            interSec.S12(countS1)=interSec.S10{i,j}(2);
            interSec.S21(countS1)=interSec.S20{i,j}(1);
            interSec.S22(countS1)=interSec.S20{i,j}(2);
        end
    end
end



%Combine consecutive intersection segments
if countS1>0
    interSec.flag0=1;
    %sort the intervals
    [interSec.S11, ind]=sort(interSec.S11);
    interSec.S12=interSec.S12(ind);
    [interSec.S21, ind]=sort(interSec.S21);
    interSec.S22=interSec.S22(ind);
    S11_merged=[];
    S12_merged=[];
    S21_merged=[];
    S22_merged=[];
    %On path 1
    starti=1;
    endi=1;
    k=1;
    S11_merged(k)=interSec.S11(starti);
    S12_merged(k)=interSec.S12(endi);
    while starti<countS1 && endi<countS1
        j=1;
        while interSec.S11(starti+j)<=interSec.S12(endi)
            endi=endi+1;
            j=j+1;
            if endi>=countS1
                break;
            end
        end
        S11_merged(k)=interSec.S11(starti);
        S12_merged(k)=interSec.S12(endi);
        starti=starti+j;
        endi=starti;
        k=k+1;
    end
    interSec.S11_merged=S11_merged;
    interSec.S12_merged=S12_merged;
    
    
    %On path 2
    starti=1;
    endi=1;
    k=1;
    S21_merged(k)=interSec.S21(starti);
    S22_merged(k)=interSec.S22(endi);
    while starti<countS1 && endi<countS1
        j=1;
        while interSec.S21(starti+j)<=interSec.S22(endi)
            endi=endi+1;
            j=j+1;
            if endi>=countS1
                break;
            end
        end
        S21_merged(k)=interSec.S21(starti);
        S22_merged(k)=interSec.S22(endi);
        starti=starti+j;
        endi=starti;
        k=k+1;
    end
    interSec.S21_merged=S21_merged;
    interSec.S22_merged=S22_merged;
end

%Find the time intervals with intersections between the two paths for the
%given velocity profiles

% v=pathVel.v;
% u_maxD=pathVel.u_maxD;
% v_maxD=pathVel.v_maxD;
% v_maxDC=pathVel.v_maxDC;
% s_bar1=pathVel.s_bar1;
% s_bar2=pathVel.s_bar2;
% for i=1:length(S11_merged)
%     Ti(i)=findTimeOnPath(S,Path,pathVel)
% end
end






