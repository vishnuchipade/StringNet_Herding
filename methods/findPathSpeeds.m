function [pathVel] = findPathSpeeds(path, v_maxDC, u_maxD, v_maxD)
% ------------------------------
% Author: Vishnu S. Chipade
% Email: vishnuc@umich.edu
% ------------------------------

% This function finds the optimal terminal speeds on each segment of the
% path and also gives times at which the agent would be at that
% intersection (quadratic drag in the double integrator dynamics)

%rI=Path.rV(:,1);
%rF=Path.rV(:,Path.NS);
global rho_safe C_d


lambda0=sqrt(C_d^2*rho_safe^2+4);
kappa1=rho_safe*u_maxD*lambda0;
kappa2=C_d*rho_safe^2*u_maxD;

lambda1=sqrt(rho_safe*(lambda0-rho_safe*C_d));
lambda2=sqrt(rho_safe*(lambda0+rho_safe*C_d));
lambda3=lambda0/rho_safe*sqrt(u_maxD/2);

S=path.S;
P=path.P;
%Number of segments
NS=length(S)-1;
if NS~=1  %If NS~=1 then there is at least one circular arc segment in the path and at least 3 total segments
    %From the beginning
    i=1;
    v(i)=0;
    while v(i)<v_maxDC
        if mod(i,2)==1
            if P(i)>1/(2*C_d)*log((u_maxD-C_d*v(i)^2)/(u_maxD-C_d*v_maxDC^2))
                v(i+1)=v_maxDC;
                lambda=(u_maxD+C_d*v(i+1)^2)/(u_maxD-C_d*v(i)^2)*exp(2*C_d*P(i));
                v_bar(i)=sqrt((lambda-1)/(lambda+1)*u_maxD/C_d);
                 if abs(v_bar(i)-v_maxD)<1e-10
                    v_bar(i)=(1-1e-10)*v_maxD;
                 end
                s_bar1(i)=1/(2*C_d)*log((u_maxD-C_d*v(i)^2)/(u_maxD-C_d*v_bar(i)^2));%abs(v_bar(i)^2-v(i)^2)/(2*u_maxD);
                s_bar2(i)=P(i)-1/(2*C_d)*log((u_maxD+C_d*v_bar(i)^2)/(u_maxD+C_d*v(i+1)^2));%%P(i)-abs(v_bar(i)^2-v_maxDC^2)/(2*u_maxD);
            else
                v(i+1)=sqrt((u_maxD-exp(-2*C_d*P(i))*(u_maxD-C_d*v(i)^2))/C_d);
                v_bar(i)=v(i+1);
                s_bar1(i)=P(i);
                s_bar2(i)=s_bar1(i);
            end
        else
            %             if S(i)>R/2*(atanh(ve^2/vmc^2)-atanh(v(i-1)^2/vmc^2))
            %                 v(i)=v_maxDC;
            %             else
            v(i+1)=sqrt(0.5*(kappa1*tanh(P(i)*lambda0/rho_safe+atanh((kappa2+2*v(i)^2)/kappa1))-kappa2));
            v_bar(i)=v(i+1);
            s_bar1(i)=P(i);
            s_bar2(i)=s_bar1(i);
            %             end
        end
        i=i+1;
    end
    i1=i;
    
    %From the end
    i=NS+1;
    v(i)=0;
    while v(i)<v_maxDC
        if mod(i,2)==0
            if P(i-1)>1/(2*C_d)*log((u_maxD+C_d*v_maxDC^2)/(u_maxD+C_d*v(i)^2)) %abs(v_maxDC^2-v(i)^2)/(2*u_maxD)
                v(i-1)=v_maxDC;
                lambda=(u_maxD+C_d*v(i)^2)/(u_maxD-C_d*v(i-1)^2)*exp(2*C_d*P(i-1));
                v_bar(i-1)=sqrt((lambda-1)/(lambda+1)*u_maxD/C_d);
                %v_bar(i-1)=0.5*sqrt(2*v(i-1)^2+2*v(i)^2+4*P(i-1)*u_maxD);
                 if abs(v_bar(i-1)-v_maxD)<1e-10
                    v_bar(i-1)=(1-1e-10)*v_maxD;
                 end
                s_bar1(i-1)=1/(2*C_d)*log((u_maxD-C_d*v(i-1)^2)/(u_maxD-C_d*v_bar(i-1)^2));%abs(v_bar(i-1)^2-v_maxDC^2)/(2*u_maxD);
                s_bar2(i-1)=P(i-1)-1/(2*C_d)*log((u_maxD+C_d*v_bar(i-1)^2)/(u_maxD+C_d*v(i)^2));%P(i-1)-abs(v_bar(i-1)^2-v(i)^2)/(2*u_maxD);
            else
                v(i-1)=sqrt((-u_maxD+exp(-2*C_d*P(i))*(u_maxD+C_d*v(i)^2))/C_d);
                v_bar(i-1)=v(i-1);
                s_bar1(i-1)=0;
                s_bar2(i-1)=s_bar1(i-1);
            end
        else
            %             lambda0=sqrt(C_d^2*rho_safe^2+4);
            %             kappa1=rho_safe*u_maxD*lambda0;
            %             kappa2=C_d*rho_safe^2*u_maxD;
            %v(i-1)=v_maxDC*sqrt(tanh(2*P(i-1)/rho_safe+atanh(v(i)^2/v_maxDC^2)));
            v(i-1)=sqrt(0.5*(kappa2-kappa1*tanh(-P(i)*lambda0/rho_safe+atanh((kappa2-2*v(i)^2)/kappa1))));
            v_bar(i-1)=v(i-1);
            s_bar1(i-1)=P(i);
            s_bar2(i-1)=s_bar1(i-1);
            %             end
        end
        i=i-1;
    end
    i2=i;
    
    %for rest of the nodes
    v(i1+1:i2-1)=ones(i2-i1-1,1)*v_maxDC;
    for i=i1:i2-1
        if mod(i,2)==1 %straight line segments
            lambda=(u_maxD-C_d*v(i+1)^2)/(u_maxD-C_d*v(i)^2)*exp(2*C_d*P(i));
            v_bar(i)=sqrt((lambda-1)/(lambda+1)*u_maxD/C_d);
            if abs(v_bar(i)-v_maxD)<1e-10
                v_bar(i)=(1-1e-10)*v_maxD;
            end
            s_bar1(i)=1/(2*C_d)*log((u_maxD-C_d*v_maxDC^2)/(u_maxD-C_d*v_bar(i)^2));%abs(v_bar(i)^2-v(i)^2)/(2*u_maxD);
            s_bar2(i)=P(i)-1/(2*C_d)*log((u_maxD+C_d*v_bar(i)^2)/(u_maxD+C_d*v_maxDC^2));
        else %circular segments
            v_bar(i)=v_maxDC;
            s_bar1(i)=0;
            s_bar2(i)=P(i);
        end
    end
    
else %Only one segment i.e. no circular collision segment (i.e no obstacle in the way)
    v=[0,0];
    i=1;
    lambda=(u_maxD+C_d*v(i+1)^2)/(u_maxD-C_d*v(i)^2)*exp(2*C_d*P(i));
    v_bar(i)=sqrt((lambda-1)/(lambda+1)*u_maxD/C_d);
    if abs(v_bar(i)-v_maxD)<1e-10
        v_bar(i)=(1-1e-10)*v_maxD;
    end
    s_bar1(i)=1/(2*C_d)*log((u_maxD)/(u_maxD-C_d*v_bar(i)^2));
    s_bar2(i)=P(i)-1/(2*C_d)*log((u_maxD+C_d*v_bar(i)^2)/(u_maxD));
end

%Return speed profiles
pathVel.v=v;
pathVel.v_bar=v_bar;
pathVel.s_bar1=s_bar1;
pathVel.s_bar2=s_bar2;

%Find the terminal times at each segment
T=zeros(1,NS+1);
T_bar1=zeros(1,NS);
T_bar2=zeros(1,NS);

T(1)=0; %Initial position

for i=1:NS
    if mod(i,2)==1 %for straight segments
%         if abs(v_bar(i)-v_maxD)>1e-10
%             T(i+1)=T(i)+1/sqrt(u_maxD*C_d)*(atanh(v_bar(i)/v_maxD)+atan((v_bar(i)/v_maxD))-atanh(v(i)/v_maxD)-atan((v(i+1)/v_maxD)));
%         else
%             T(i+1)=T(i)+1/sqrt(u_maxD*C_d)*(atan((v_bar(i)/v_maxD))-atan((v(i+1)/v_maxD)));
%         end
        
        T_bar1(i)=T(i)+1/sqrt(u_maxD*C_d)*(atanh(v_bar(i)/v_maxD)-atanh(v(i)/v_maxD));
        T_bar2(i)=T_bar1(i)+(s_bar2(i)-s_bar1(i))/v_bar(i);
        T(i+1)= T_bar2(i) + 1/sqrt(u_maxD*C_d)*(atan((v_bar(i)/v_maxD))-atan((v(i+1)/v_maxD)));
        %         if P(i)<0.5*((v_maxD^2-v(i)^2)+(v_maxD^2-v(i+1)^2))/u_maxD
        %             T(i+1)=T(i) + (2*v_bar(i)-v(i)-v(i+1))/u_maxD;
        %         else
        %             T(i+1)=T(i) + P(i)/v_maxD-0.5*((v_maxD^2-v(i)^2)+(v_maxD^2-v(i+1)^2))/(u_maxD*v_maxD);
        %         end
    else %For circular arc segments
        if v(i)<v_maxDC && v(i+1)<v_maxDC && v(i)<v(i+1)
            
            %             kappa0=rho_safe/2/lambda0;
            %             kappa=P(i)+2*kappa0*(atanh((kappa2-2*v(i)^2)/kappa1)-atanh((kappa2-2*v(i+1)^2)/kappa1));
            %             ek=exp(kappa/kappa0);
            
            T(i+1)=T(i) +1/lambda3*(atan(sqrt(2/u_maxD)*v(i+1)/lambda1)/lambda1+atanh(sqrt(2/u_maxD)*v(i+1)/lambda2)/lambda2)...
                -1/lambda3*(atan(sqrt(2/u_maxD)*v(i)/lambda1)/lambda1+atanh(sqrt(2/u_maxD)*v(i)/lambda2)/lambda2) ;
            T_bar1(i)=T(i+1);
            T_bar2(i)=T_bar1(i);
        elseif v(i)<v_maxDC && v(i+1)<v_maxDC && v(i)>v(i+1)
            T(i+1)=T(i) +1/lambda3*(atan(sqrt(2/u_maxD)*v(i+1)/lambda2)/lambda2+atanh(sqrt(2/u_maxD)*v(i+1)/lambda1)/lambda1)...
                -1/lambda3*(atan(sqrt(2/u_maxD)*v(i)/lambda2)/lambda2+atanh(sqrt(2/u_maxD)*v(i)/lambda1)/lambda1) ;
            T_bar1(i)=T(i);
            T_bar2(i)=T_bar1(i);
        else
            
            T(i+1)=T(i) + P(i)/v_maxDC;
            T_bar1(i)=T(i);
            T_bar2(i)=T(i+1);
        end
    end
    
end

%return time profile as well
pathVel.T=T;
pathVel.T_bar1=T_bar1;
pathVel.T_bar2=T_bar2;
pathVel.u_maxD=u_maxD;
pathVel.v_maxD=v_maxD;
pathVel.v_maxDC=v_maxDC;
end