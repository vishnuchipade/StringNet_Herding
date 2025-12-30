function [interSec]=pathIntersections(Path1,Path2)
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
        
        interSec0=interSecLineLine(r1,r2,mL1,cL1,theta1,drx,dry,L1,r21,r22,mL2,cL2,theta2,drx2,dry2,L2,dtheta);
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
    
    %with circular segments of the second path
    for j=2:2:NS2
        Sj0=Path2.S(j);
        r21=rV2(:,j);
        r22=rV2(:,j+1);
        interSec0=interSecLineCirc(r1,r2,mL1,cL1,L1,dr_hat,drx,dry,r21,r22,rVC2(:,j),segType2(j));
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
        interSec0=interSecLineCirc(r1,r2,mL1,cL1,L1,dr_hat,drx,dry,r21,r22,rVC1(:,i),segType1(i));
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
        Sj0=Path2.S(j);
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
                ds1=rho_safe*acos(1-norm(interSec.Pos1{i,j}(:,1)-interSec.Pos1{i,j}(:,2))^2/(2*rho_safe^2));
                interSec.S10{i,j}(1)=Si0+rho_safe*acos(1-norm(interSec.Pos1{i,j}(:,1)-r1)^2/(2*rho_safe^2));
                interSec.S10{i,j}(2)=interSec.S10{i,j}(1)+ds1;
                interSec.Pbar1=interSec.Pbar1+ds1;
                
                interSec.Pos2{i,j}(:,1)=rVC1(:,i)+rho_safe*[cos(ang11+angfi(3)), sin(ang11+angfi(3))]';  %middle angles will give the end points of the common segment
                interSec.Pos2{i,j}(:,2)=rVC1(:,i)+rho_safe*[cos(ang11+angfi(4)), sin(ang11+angfi(4))]';
                ds2=rho_safe*acos(1-norm(interSec.Pos2{i,j}(:,1)-interSec.Pos2{i,j}(:,2))^2/(2*rho_safe^2));
                interSec.S20{i,j}(1)=Sj0+rho_safe*acos(1-norm(interSec.Pos2{i,j}(:,1)-r21)^2/(2*rho_safe^2));
                interSec.S20{i,j}(2)=interSec.S20{i,j}(1)+ds2;
                interSec.Pbar2=interSec.Pbar2+ds2;
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
    
    Si11=[];
    Si12=[];
    Si21=[];
    Si22=[];
    %On path 1
    starti=1;
    endi=1;
    k=1;
    Si11(k)=interSec.S11(starti);
    Si12(k)=interSec.S12(endi);
    while starti<countS1 && endi<countS1
        j=1;
        while interSec.S11(starti+j)-interSec.S12(endi)<=1e-10  %(use a small positive number instead of zero)
            endi=endi+1;
            j=j+1;
            if endi>=countS1
                break;
            end
        end
        Si11(k)=interSec.S11(starti);
        Si12(k)=interSec.S12(endi);
        starti=starti+j;
        endi=starti;
        k=k+1;
    end
    interSec.Si11=Si11;
    interSec.Si12=Si12;
    
    
    %On path 2
    starti=1;
    endi=1;
    k=1;
    Si21(k)=interSec.S21(starti);
    Si22(k)=interSec.S22(endi);
    while starti<countS1 && endi<countS1
        j=1;
        while interSec.S21(starti+j)-interSec.S22(endi)<=1e-10 %(use a small positive number instead of zero)
            endi=endi+1;
            j=j+1;
            if endi>=countS1
                break;
            end
        end
        Si21(k)=interSec.S21(starti);
        Si22(k)=interSec.S22(endi);
        starti=starti+j;
        endi=starti;
        k=k+1;
    end
    interSec.Si21=Si21;
    interSec.Si22=Si22;
end

%Find the time intervals with intersections between the two paths for the
%given velocity profiles

% v=pathVel.v;
% u_maxD=pathVel.u_maxD;
% v_maxD=pathVel.v_maxD;
% v_maxDC=pathVel.v_maxDC;
% s_bar1=pathVel.s_bar1;
% s_bar2=pathVel.s_bar2;
% for i=1:length(Si11)
%     Ti(i)=findTimeOnPath(S,Path,pathVel)
% end
end






