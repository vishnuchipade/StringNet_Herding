function [Psi, Psi_dot, Psi_ddot] = formationRefTraj(rA,rA_dot,rA_ddot,rO,rS)
% ------------------------------
% Author: Vishnu S. Chipade
% Email: vishnuc@umich.edu
% ------------------------------

%This function gives direction using vector fields based on super-elliptic contours around rectangular
%obstclaes
%rF=[xF,yF]';  rO=[xO1,xO2,aO1,bO1;xO2,xO2,aO2,bO2;....]';
%rS=[xS,yS]';
%rho_Fmax is the radius of the formation
global E_m_O E_u_O E_bar_O A B C D aO bO nO sigmaProd_dot_fun sigmaProd_ddot_fun d2Beta_barF_dBeta2_fun
NO=size(rO,2);
NA=size(rA,2);
for i=1:NA
    rF=rA(:,i);
    rF_dot=rA_dot(:,i);
    rF_ddot=rA_ddot(:,i);
%Based on potential fields (rectangular obstacles)
%F_S=zeros(2,1);
R_SF=norm(rS-rF);
R_SF_dot=-(rS-rF)'*rF_dot/R_SF;
R_SF_ddot=(-(rS-rF)'*rF_dot+norm(rF_dot)^2-R_SF_dot^2)/R_SF;
F_O=zeros(2,1);
%R2=sqrt((rS(1)-rF(1))^2+(rS(2)-rF(2))^2);
sigmaProd=1;
Sigma=zeros(NO,1);
Sigma_dot=Sigma;
Sigma_ddot=Sigma;
F_Ok=zeros(2,NO);
F_Ok_dot=F_Ok;
F_Ok_ddot=F_Ok;
for k=1:NO
    n=nO(k);
    a=aO(k);
    b=bO(k);
    R_OF=norm(rF(:)-rO(1:2,k));
    R_FS=norm(rS-rF);
    R_OF_dot=(rF(:)-rO(1:2,k))'*rF_dot/R_OF;
    R_FS_dot=-(rS(:)-rF(1:2))'*rF_dot/R_FS;
    E_Ok=(abs((rF(1)-rO(1,k))/a))^(2*n)+(abs((rF(2)-rO(2,k))/b))^(2*n)-1;
    E_Ok_dot=2*n*[(abs((rF(1)-rO(1,k))/a))^(2*n)/(rF(1)-rO(1,k)), (abs((rF(2)-rO(2,k))/b))^(2*n)/(rF(2)-rO(2,k))]*rF_dot;
    E_Ok_ddot=2*n*(n-1)*[(abs((rF(1)-rO(1,k))/a))^(2*n)/(rF(1)-rO(1,k))^2,  (abs((rF(2)-rO(2,k))/b))^(2*n)/(rF(2)-rO(2,k))^2]*rF_dot.^2 ...
             + 2*n*[(abs((rF(1)-rO(1,k))/a))^(2*n)/(rF(1)-rO(1,k)),  (abs((rF(2)-rO(2,k))/b))^(2*n)/(rF(2)-rO(2,k))]*rF_ddot;
    if E_Ok>E_m_O(k) && E_Ok<E_bar_O(k)
        sigma=1;
        sigma_dot=0;    
        sigma_ddot=0;
    elseif E_Ok>E_bar_O(k) && E_Ok<E_u_O(k)
        sigma=A(k)*E_Ok^3+B(k)*E_Ok^2+C(k)*E_Ok+D(k);
        sigma_dot=(3*A(k)*E_Ok^2+2*B(k)*E_Ok+C(k))*E_Ok_dot;%h=1;
        sigma_ddot=(6*A(k)*E_Ok+2*B(k))*(E_Ok_dot)^2+(3*A(k)*E_Ok^2+2*B(k)*E_Ok+C(k))*E_Ok_ddot;
    else
        sigma=0;
        sigma_dot=0;
        sigma_ddot=0;
    end
    Sigma(k)=sigma;
    Sigma_dot(k)=sigma_dot;
    Sigma_ddot(k)=sigma_ddot;
    %E_O=E_Ok;
    betaOS=atan2(rS(2)-rO(2,k),rS(1)-rO(1,k));  %angle between the safe area and the obstacle
    betaOF=atan2(rF(2)-rO(2,k),rF(1)-rO(1,k));   %angle between the location and the obstacle
    betaOF_dot=[-(rF(2)-rO(2,k)), (rF(1)-rO(1,k))]*rF_dot/R_OF^2;
    betaOF_ddot=([-(rF(2)-rO(2,k)), (rF(1)-rO(1,k))]*rF_ddot-2*R_OF*R_OF_dot*betaOF_dot)/R_OF^2;
    beta_barOS=atan2(b^(2*n)*sign(rS(1)-rO(1,k))*(abs(rS(1)-rO(1,k)))^(2*n-1),-a^(2*n)*sign(rS(2)-rO(2,k))*(abs(rS(2)-rO(2,k)))^(2*n-1));  %angle of the tangent at the safe area
    betaSF=atan2(rF(2)-rS(2),rF(1)-rS(1));
    betaSF_dot=[-(rF(2)-rS(2)), (rF(1)-rS(1))]*(rF_dot)/R_FS^2;
    betaSF_ddot=([-(rF(2)-rS(2)), (rF(1)-rS(1))]*rF_ddot-2*R_FS*R_FS_dot*betaOF_dot)/R_FS^2;
            
    if E_Ok>E_m_O(k) 
        beta_bar=atan2(b^(2*n)*sign(rF(1)-rO(1,k))*(abs(rF(1)-rO(1,k)))^(2*n-1),-a^(2*n)*sign(rF(2)-rO(2,k))*(abs(rF(2)-rO(2,k)))^(2*n-1));  %angle of the tangent at the given location
        if beta_bar<0   % get betaOS between [0,2pi]
            beta_bar=beta_bar+2*pi;
        end
        
        Delta_beta=betaOF-betaOS;
        if Delta_beta<0
            Delta_beta=Delta_beta+2*pi;
        end
        
          Delta_beta_bar=beta_barOS-betaOS;
        if Delta_beta_bar<0
            Delta_beta_bar=Delta_beta_bar+2*pi;
        end
        
         crossProd=cross([cos(betaSF);sin(betaSF);0],[cos(beta_bar);sin(beta_bar);0]);
            if crossProd(3)<0
                beta_bar1=betaSF-pi;
                beta_bar1_dot=betaSF_dot;
                beta_bar1_ddot=betaSF_ddot;
            else
                if Delta_beta>pi
                    %                 beta_bar1=beta_bar-Delta_beta_bar+Delta_beta/pi*(Delta_beta_bar-pi);
                    %                 beta_bar1_dot=(sin(2*betaOA)/(0.5*((cos(2*betaOA))^2+1)) + (Delta_beta_bar-pi)/pi)*betaOA_dot;
                    beta_bar1=beta_bar;
                else
                    %                 beta_bar1=beta_bar-Delta_beta_bar/pi*(Delta_beta-pi);
                    %                 beta_bar1_dot=(sin(2*betaOA)/(0.5*((cos(2*betaOA))^2+1))-(Delta_beta_bar)/pi)*betaOA_dot;
                    beta_bar1=beta_bar-pi;
                end
                 s0=num2cell([betaOF,n,a,b]');
                beta_bar1_dot=(b^(2*n)*(2*n-1)*sign(cos(betaOF))*abs(cos(betaOF))^(2*n-2))/(a^(2*n)*sign(sin(betaOF))*abs(sin(betaOF))^(2*n))*betaOF_dot;
                beta_bar1_ddot=(b^(2*n)*(2*n-1)*sign(cos(betaOF))*abs(cos(betaOF))^(2*n-2))/(a^(2*n)*sign(sin(betaOF))*abs(sin(betaOF))^(2*n))*betaOF_ddot...
                   +d2Beta_barF_dBeta2_fun(s0{:})*betaOF_dot^2;
            end
        
        F_Ok(:,k)=[cos(beta_bar1); sin(beta_bar1)]; 
        F_Ok_dot(:,k)=beta_bar1_dot*[-sin(beta_bar1);cos(beta_bar1)];
        F_Ok_ddot(:,k)=beta_bar1_ddot*[-sin(beta_bar1);cos(beta_bar1)]-beta_bar1_dot^2*F_Ok(:,k);
        F_O=F_O+sigma*F_Ok(:,k);
        sigmaProd=sigmaProd*(1-sigma);
    end
    
end
F_S=(rS-rF)/R_SF;
F=F_O+sigmaProd*F_S;
F_S_dot=(-rF_dot-F_S*R_SF_dot)/R_SF;
F_S_ddot=(-rF_ddot-R_SF_ddot*F_S-2*R_SF_dot*F_S_dot)/R_SF;
%to find F_dot
F_O_dot=0;
F_O_ddot=0;
sigmaProd_dot=0;
sigmaProd_ddot=0;
for k=1:NO
    F_O_dot=F_O_dot+Sigma(k)*F_Ok_dot(:,k)+Sigma_dot(k)*F_Ok(:,k);
    F_O_ddot=F_O_ddot+Sigma(k)*F_Ok_ddot(:,k)+2*Sigma_dot(k)*F_Ok_dot(:,k)+Sigma_ddot(k)*F_Ok(:,k);
    temp3=1;
    for kk=1:NO
        if kk~=k
        temp3=temp3*(1-Sigma(kk));
        else
            temp3=-Sigma_dot(kk)*temp3;
        end
    end
   sigmaProd_dot=sigmaProd_dot+temp3;
end
s1=num2cell([Sigma',Sigma_dot']);
s2=num2cell([Sigma',Sigma_dot',Sigma_ddot']);
sigmaProd_dot=sigmaProd_dot_fun(s1{:});
sigmaProd_ddot=sigmaProd_ddot_fun(s2{:});
F_dot=sigmaProd_dot*F_S+sigmaProd*F_S_dot+F_O_dot;
F_ddot=sigmaProd_ddot*F_S+2*sigmaProd_dot*F_S_dot+sigmaProd*F_S_ddot+F_O_ddot;
psi=atan2(F(2),F(1));
norm_F=norm(F);
psi_dot=[-F(2),F(1)]*F_dot/norm_F^2;
psi_ddot=([-F(2),F(1)]*F_ddot -2*sum(F.*F_dot)*psi_dot)/norm_F;
Psi(i,1)=psi;
Psi_dot(i,1)=psi_dot;
Psi_ddot(i,1)=psi_ddot;
end
end


%     if betaOS<0   % get betaOS between [0,2pi]
%         betaOS=betaOS+2*pi;
%     end
%     if betaOF<0   % get betaOF between [0,2pi]
%         betaOF=betaOF+2*pi;
%     end
%     if beta_barOS<0   % get beta_barOS between [0,2pi]
%         beta_barOS=beta_barOS+2*pi;
%     end


% if (1)
%         if betaOS<pi
%             if betaOF<=betaOS
%                 beta_bar1=beta_bar-((-betaOS+pi+betaOF)/pi)*(beta_barOS-betaOS);
%             elseif betaOF>betaOS+pi
%                 % beta_bar1=betaOS+((betaOS-beta0+pi)/pi)*(beta0-betaOS);
%                 beta_bar1=beta_bar+((betaOS+pi-betaOF)/pi)*(beta_barOS-betaOS);
%             else
%                 % beta_bar1=betaOS+((beta0-betaOS)/pi)*(betaOF-betaOS);
%                 
%                 beta_bar1=beta_bar-(beta_barOS-betaOS)+(beta_barOS-betaOS-pi)*(betaOF-betaOS)/pi;
%             end
%         elseif betaOS>=pi && betaOS<1.5*pi
%             if betaOF<=betaOS && betaOF>betaOS-pi
%                 beta_bar1=beta_bar+((betaOS-betaOF-pi)/pi)*(beta_barOS-betaOS);
%             elseif betaOF>betaOS
%                 beta_bar1=beta_bar-(beta_barOS-betaOS)+(beta_barOS-betaOS-pi)*(betaOF-betaOS)/pi;
%             else
%                 beta_bar1=beta_bar-(beta_barOS-betaOS)+(beta_barOS-betaOS-pi)*(betaOF+2*pi-betaOS)/pi;
%             end
%         else
%             if beta_barOS<pi/2
%                 beta_barOS=beta_barOS+2*pi;
%             end
%             if betaOF<=betaOS && betaOF>betaOS-pi
%                 beta_bar1=beta_bar+((betaOS-betaOF-pi)/pi)*(beta_barOS-betaOS);
%             elseif betaOF>betaOS
%                 beta_bar1=beta_bar-(beta_barOS-betaOS)+(beta_barOS-betaOS-pi)*(betaOF-betaOS)/pi;
%             else
%                 beta_bar1=beta_bar-(beta_barOS-betaOS)+(beta_barOS-betaOS-pi)*(betaOF+2*pi-betaOS)/pi;
%             end
%         end
%         end
