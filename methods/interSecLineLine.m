function interSec=interSecLineLine(r1,r2,mL1,cL1,theta1,drx,dry,L1,r21,r22,mL2,cL2,theta2,drx2,dry2,L2,dtheta)
% ------------------------------
% Author: Vishnu S. Chipade
% Email: vishnuc@umich.edu
% ------------------------------


global rho_D

interSec.Flag=0;
if mL1==mL2
    %&& (r1(1)==r21(1) && r1(2)==r21(2)) || (r1(1)==r22(1) && r1(2)==r22(2))
    %find the distance between the two || lines
    d=abs(cL1-cL2)/sqrt(1+mL1^2);
    if d<2*rho_D  %intersection of the rectangular area swept by the footprint of the defender
        interSec.Flag=1;   
        if L1<L2
            interSec.Pos1(:,1)=r1;  %first end point of the intersection segment on path 1
            interSec.Pos1(:,2)=r2;   %second end point of the intersection segment on path 1
            interSec.S10(1)=0;
            interSec.S10(2)=L1;
            interSec.Pbar1=L1;
            
            interSec.Pos2(1,1)=(r1(1)+mL1*r1(2)-mL1*cL2)/(1+mL1^2);  %x-coord of first end point of the intersection segment on path 2
            interSec.Pos2(1,2)=(r2(1)+mL1*r2(2)-mL1*cL2)/(1+mL1^2);  %x-coord of second end point of the intersection segment on path 2
            if mL2<1e16
                interSec.Pos2(2,1)=mL2*interSec.Pos2(1,1)+cL2;
                interSec.Pos2(2,2)=mL2*interSec.Pos2(1,2)+cL2;
            else
                interSec.Pos2(2,1)=r1(2);
                interSec.Pos2(2,2)=r2(2);
            end
            interSec.S20(1)=norm(r21-interSec.Pos2(:,1));
            interSec.S20(2)=interSec.S20(1)+L1;
            interSec.Pbar2=L1;
        else
            interSec.Pos2(:,1)=r21;  %first end point of the intersection segment on path 2
            interSec.Pos2(:,2)=r22;   %second end point of the intersection segment on path 2
            interSec.S20(1)=0;
            interSec.S20(2)=L2;
            interSec.Pbar2=L2;
            
            interSec.Pos1(1,1)=(r21(1)+mL1*r21(2)-mL1*cL2)/(1+mL1^2);  %x-coord of  first end point of the intersection segment on path 1
            interSec.Pos1(1,2)=(r22(1)+mL1*r22(2)-mL1*cL2)/(1+mL1^2);  %x-coord of  second end point of the intersection segment on path 1
            if mL1<1e16
                interSec.Pos1(2,1)=mL1*interSec.Pos1(1,1)+cL1;
                interSec.Pos1(2,2)=mL1*interSec.Pos1(1,2)+cL1;
            else
                interSec.Pos1(2,1)=r21(2);
                interSec.Pos1(2,2)=r22(2);
            end
            interSec.S10(1)=norm(r1-interSec.Pos1(:,1));
            interSec.S10(2)=interSec.S10(1)+L2;
            interSec.Pbar1=L2;
        end
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
    %interSec.Flag=((0<= lambda11) && (1>= lambda11)) || ((0<= lambda12) && (1>= lambda12)) || ...
    %((0<= lambda21) && (1>= lambda21)) ||  ((0<= lambda22) && (1>= lambda22));
    
    interSec.Flag=0;
    %             if (0<= lambda1) && (0<= lambda2)   %both path's end part are possibly in conflict
    %                 if (lambda11<=1 || lambda12<=1) && (lambda21<=1 || lambda22<=1)
    %                     interSec.Flag=1;
    %                 end
    %             elseif (0<= lambda1) && (0> lambda2)
    %                 if (lambda11<=1 || lambda12<=1) && (lambda21>0 || lambda22>=0)
    %                     interSec.Flag=1;
    %                 end
    %             elseif (0> lambda1) && (0<= lambda2)
    %                 if (lambda11>0 || lambda12>=0) && (lambda21<=1 || lambda22<=1)
    %                     interSec.Flag=1;
    %                 end
    %             else %(0> lambda1) && (0> lambda2)
    %                 if (lambda11>0 || lambda12>=0) && (lambda21>0 || lambda22>=0)
    %                     interSec.Flag=1;
    %                 end
    %             end
    if lambda2>=1
        if lambda1>0 && lambda1<1
            % interSec.Flag=1;   %1 if point intersection
            if norm(r22-[x_int,y_int]')<2*rho_D/abs(sin(dtheta))
                interSec.Flag=1;
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
                        if min(lambda00)<=1
                        lam21=max(0,min(lambda00));
                        else
                            lam21=1;
                             interSec.Flag=1;
                        end
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
                         if min(lambda00)<=1
                        lam21=max(0,min(lambda00));
                        else
                            lam21=1;
                             interSec.Flag=1;
                        end
                    else
                        lam11=lambda11;
                        lam21=max(0,lambda21);
                    end
                end
                
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
                
                interSec.Flag=1;
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
                
                interSec.Flag=1;
            end
        end
    elseif lambda2>0 && lambda2<1
        if lambda1>0 && lambda1<1
            %interSec.Flag=1;   %1 if point intersection
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
            interSec.Flag=1;
        elseif lambda1<=0
            if norm(r1-[x_int,y_int]')<2*rho_D/abs(sin(dtheta))
                 interSec.Flag=1;
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
                        if max(lambda00)>=0
                        lam12=min(1,max(lambda00));
                        else
                            lam12=0;
                             interSec.Flag=0;
                        end
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
                         if max(lambda00)>=0
                        lam12=min(1,max(lambda00));
                        else
                            lam12=0;
                             interSec.Flag=0;
                        end
                    else
                        lam21=lambda21;
                        lam12=min(1,lambda12);
                    end
                end
               
            end
        elseif lambda1>=1
            if norm(r2-[x_int,y_int]')<2*rho_D/abs(sin(dtheta))
                interSec.Flag=1;
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
                        
                        if min(lambda00)<=1
                        lam11=max(0,min(lambda00));
                        else
                            lam11=1;
                            interSec.Flag=0;
                        end
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
                        if min(lambda00)<=1
                        lam11=max(0,min(lambda00));
                        else
                            lam11=1;
                            interSec.Flag=0;
                        end
                    else
                        lam21=lambda21;
                        lam11=max(0,lambda11);
                    end
                end
                
            end
        end
    else %lambda2<0
        if lambda1>0 && lambda1<1
            %interSec.Flag=1;   %1 if point intersection
            if norm(r21-[x_int,y_int]')<2*rho_D/abs(sin(dtheta))
                 interSec.Flag=1;
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
                        if max(lambda00)>=0
                        lam22=min(1,max(lambda00));
                        else
                            lam22=0;
                            interSec.Flag=0;
                        end
                        
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
                         if max(lambda00)>=0
                        lam22=min(1,max(lambda00));
                        else
                            lam22=0;
                            interSec.Flag=0;
                        end
                    else
                        lam11=lambda11;
                        lam22=min(1,lambda22);
                    end
                end
               
            end
        elseif lambda1<=0
            if norm(r21-r1)<2*rho_D
                lam21=0;
                lam11=0;
                if (r22-[x_int,y_int]')'*(r1-[x_int,y_int]')>0
                    lam12=min(1,lambda12);
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
                interSec.Flag=1;
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
                interSec.Flag=1;
            end
        end
    end
    
    if interSec.Flag==1
        interSec.Pos1(:,1)=[r1(1)+lam11*drx, r1(2)+lam11*dry]';
        interSec.Pos1(:,2)=[r1(1)+lam12*drx, r1(2)+lam12*dry]';
        interSec.Pbar1=(lam12-lam11)*L1;
        interSec.S10(1)=lam11*L1;
        interSec.S10(2)=lam12*L1;
        
        interSec.Pos2(:,1)=[r21(1)+lam21*drx2, r21(2)+lam21*dry2]';
        interSec.Pos2(:,2)=[r21(1)+lam22*drx2, r21(2)+lam22*dry2]';
        interSec.Pbar2=(lam22-lam21)*L2;
        interSec.S20(1)=lam21*L2;
        interSec.S20(2)=lam22*L2;
    end
    
end

end