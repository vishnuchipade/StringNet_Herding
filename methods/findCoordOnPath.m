function [rp, thetap]=findCoordOnPath(S,path)
% ------------------------------
% Author: Vishnu S. Chipade
% Email: vishnuc@umich.edu
% ------------------------------


global rho_safe
rp=zeros(2,length(S));
for i=1:length(S)
    %S=Si11;
    ind0=find(path.S<S(i));
    if ~isempty(ind0)
        segId=ind0(end);
    else
        segId=1;
    end    
    if segId>=length(path.S)
        1
    end
    S0=path.S(segId);
    S1=path.S(segId+1);
    dS=S(i)-S0;
    r11=path.rV(:,segId);
    r12=path.rV(:,segId+1);
    if mod(segId,2)==1  %if segment is straight line
        lambda=dS/(S1-S0);
        rp(:,i)=(1-lambda)*r11+lambda*r12;
        thetap(:,i)=atan2(r12(2)-r11(2),r12(1)-r11(1));
    else %if the segment is circular arc
        rVC=path.rVC(:,segId);
        ang1=atan2(r11(2)-rVC(2),r11(1)-rVC(1));
        if ang1<0
            ang1=ang1+2*pi;
        end
        rotMat=[cos(ang1), sin(ang1);  -sin(ang1), cos(ang1)];
        ang2=atan2(rotMat(2,:)*r12,rotMat(1,:)*r12);
        ang=dS/(S1-S0)*(ang2-ang1);
        rp(:,i)=rVC+rho_safe*[cos(ang1+ang), sin(ang1+ang)]';
        if ang<0
            ang=ang+2*pi;
        end
        if path.segType(segId)==1  %anticlockwsie segment
        thetap(:,i)=ang+pi/2;
        else
            thetap(:,i)=ang-pi/2;
        end
    end
end
end