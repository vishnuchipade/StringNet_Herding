function T=findTimeOnPath(S,Path,pathVel)
% ------------------------------
% Author: Vishnu S. Chipade
% Email: vishnuc@umich.edu
% ------------------------------

%For double integrator with quadratic drag term
global rho_safe C_d

v=pathVel.v;
u_maxD=pathVel.u_maxD;
v_maxD=pathVel.v_maxD;
v_maxDC=pathVel.v_maxDC;
s_bar1=pathVel.s_bar1;
s_bar2=pathVel.s_bar2;
T_bar1=pathVel.T_bar1;
T_bar2=pathVel.T_bar2;
v_bar=pathVel.v_bar;

lambda0=sqrt(C_d^2*rho_safe^2+4);
kappa1=rho_safe*u_maxD*lambda0;
kappa2=C_d*rho_safe^2*u_maxD;

lambda1=sqrt(rho_safe*(lambda0-rho_safe*C_d));
lambda2=sqrt(rho_safe*(lambda0+rho_safe*C_d));
lambda3=lambda0/rho_safe*sqrt(u_maxD/2);


for i=1:length(S)
    %S=Si11;
    ind0=find(Path.S<S(i));
    if ~isempty(ind0)
        segId=ind0(end);
    else
        segId=1;
    end
    S0=Path.S(segId);
    dS=S(i)-S0;
    dT_bar1=T_bar1(segId)-pathVel.T(segId);
    dT_bar2=T_bar2(segId)-pathVel.T(segId);
    if mod(segId,2)==1  %if segment is straight line
        if dS<=s_bar1(segId)
            %acceleration phase
            vtemp=(u_maxD-exp(-2*C_d*dS)*(u_maxD-C_d*v(segId)^2))/C_d;
            if vtemp<0
                v2=imag(sqrt(vtemp));  %when the term inside is very small and negative
            else
                v2=sqrt(vtemp);
            end
            dT=1/sqrt(u_maxD*C_d)*(atanh(v2/v_maxD)-atanh(v(segId)/v_maxD));
        elseif dS>s_bar1(segId) && dS<s_bar2(segId)
            dT=dT_bar1+dS/v_bar(segId);
        else            
            %acceleration plus deceleration phase
            vtemp=(-u_maxD+exp(-2*C_d*(dS-s_bar2(segId)))*(u_maxD+C_d*v_bar(segId)^2))/C_d;
             if vtemp<0
                v2=imag(sqrt(vtemp));  %when the term inside is very small and negative
            else
                v2=sqrt(vtemp);
            end
            dT=dT_bar2+ 1/sqrt(u_maxD*C_d)*(atan(v_bar(segId)/v_maxD)-atan(v2)/v_maxD);
        end
    else %if the segment is circular arc
        if dS<=s_bar1(segId)
            %acceleration phase
            v2=0.5*sqrt(kappa1*tanh(dS*lambda0/rho_safe+atanh((kappa2+2*v(segId)^2)/kappa1))-kappa2);
            dT=1/lambda3*(atan(sqrt(2/u_maxD)*v2/lambda1)/lambda1+atanh(sqrt(2/u_maxD)*v2/lambda2)/lambda2)...
                -1/lambda3*(atan(sqrt(2/u_maxD)*v(segId)/lambda1)/lambda1+atanh(sqrt(2/u_maxD)*v(segId)/lambda2)/lambda2) ;
        elseif dS>s_bar1(segId) && dS<s_bar2(segId)
            dT=dS/v_maxDC;  %terminal speeds would be same and equal to max allowed on circular arc
        else
%             %acceleration phase
%             v1_bar=v_maxDC*sqrt(tanh(2*s_bar1(segId)/rho_safe+atanh(v(segId)^2/v_maxDC^2)));
%             dT=0.5*rho_safe/v_maxDC*((atan(v1_bar/v_maxDC)+atanh(v1_bar/v_maxDC))-((atan(v(segId)/v_maxDC)+atanh(v(segId)/v_maxDC))));
%             %No zero acceleration phase due to asymptotic nature of velocity (s_bar1 and s_bar2 coincide)
            %deceleration phase
            v2=0.5*sqrt(kappa2-kappa1*tanh(-dS*lambda0/rho_safe+atanh((kappa2-2*v(segId)^2)/kappa1)));
            dT=1/lambda3*(atan(sqrt(2/u_maxD)*v2/lambda2)/lambda2+atanh(sqrt(2/u_maxD)*v2/lambda1)/lambda1)...
                -1/lambda3*(atan(sqrt(2/u_maxD)*v(segId)/lambda2)/lambda2+atanh(sqrt(2/u_maxD)*v(segId)/lambda1)/lambda1) ;
        end
    end
    T(i)=pathVel.T(segId)+dT;
end


end