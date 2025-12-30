function S=findPosOnPath(T,path,pathVel)
% ------------------------------
% Author: Vishnu S. Chipade
% Email: vishnuc@umich.edu
% ------------------------------


global rho_safe C_d

S=zeros(1,length(T));

v=pathVel.v;
v_bar=pathVel.v_bar;
u_maxD=pathVel.u_maxD;
v_maxD=pathVel.v_maxD;
v_maxDC=pathVel.v_maxDC;
s_bar1=pathVel.s_bar1;
s_bar2=pathVel.s_bar2;
v_bar=pathVel.v_bar;
T_bar1=pathVel.T_bar1;
T_bar2=pathVel.T_bar2;

lambda0=sqrt(C_d^2*rho_safe^2+4);
kappa1=rho_safe*u_maxD*lambda0;
kappa2=C_d*rho_safe^2*u_maxD;

lambda1=sqrt(rho_safe*(lambda0-rho_safe*C_d));
lambda10=sqrt(2/u_maxD)/lambda1;
lambda2=sqrt(rho_safe*(lambda0+rho_safe*C_d));
lambda20=sqrt(2/u_maxD)/lambda2;
lambda3=lambda0/rho_safe*sqrt(u_maxD/2);

for i=1:length(T)
    %S=Si11;
    ind0=find(pathVel.T<T(i));
    if ~isempty(ind0)
        segId=ind0(end);
    else
        segId=1;
    end
    if segId<length(T)
    T0=pathVel.T(segId);
    Tbar1=T_bar1(segId);
    Tbar2=T_bar2(segId);
    dT=T(i)-T0;
    if mod(segId,2)==1  %if segment is straight line
        if T(i)<=Tbar1
            %acceleration phase
            %find the speed first
            vT=sqrt(u_maxD/C_d)*tanh(sqrt(u_maxD*C_d)*dT+atanh(sqrt(C_d/u_maxD)*v(segId)));
            dS=1/(2*C_d)*log((u_maxD-C_d*v(segId)^2)/(u_maxD-C_d*vT^2));
        elseif T(i)>Tbar1 && T(i)<Tbar2
            dS=s_bar1(segId)+v_maxD*(T(i)-Tbar1);
        else
            %acceleration, zero, deceleration phase
            vT=sqrt(u_maxD/C_d)*tan(sqrt(u_maxD*C_d)*dT+atan(sqrt(C_d/u_maxD)*v_bar(segId)));
            dS=s_bar2(segId)+1/(2*C_d)*log((u_maxD+C_d*v_bar(segId)^2)/(u_maxD+C_d*vT^2));
        end
    else %if the segment is circular arc
        if T(i)<=Tbar1
            %acceleration phase
            f=@(v) 1/lambda3*(atan(lambda10*v)/lambda1 + atanh(lambda20*v)/lambda2 - atan(lambda10*v(segId))/lambda1 - atanh(lambda20*v(segId))/lambda2)-dT;
            vT=fsolve(f,v_maxDC/2);
            dS=rho_safe/lambda0*(atanh((kappa2+2*vT^2)/kappa1)-atanh((kappa2+2*v(segId)^2)/kappa1));
            %0.5*rho_safe/v_maxDC*((atan(v1/v_maxDC)+atanh(v1/v_maxDC))-((atan(v(segId)/v_maxDC)+atanh(v(segId)/v_maxDC))));
        elseif T(i)>Tbar1 && T(i)<Tbar2
            dS=dT*v_maxDC;  %terminal speeds would be same and equal to max allowed on circular arc
        else
            %acceleration phase
            f=@(v) 1/lambda3*(atan(lambda20*v)/lambda2 + atanh(lambda10*v)/lambda1 - atan(lambda20*v_bar(segId))/lambda2 - atanh(lambda10*v_bar(segId))/lambda1)-dT;
            v1=fsolve(f,v_maxDC/2);
            dS=s_bar2(segId)+rho_safe/lambda0*(atanh((kappa2-2*vT^2)/kappa1)-atanh((kappa2+2*v_bar(segId)^2)/kappa1));
        end
    end
    S(i)=path.S(segId)+dS;
    else
    S(i) = path.S(segId);
    end
end

end