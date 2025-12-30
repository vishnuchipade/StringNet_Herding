function [ uA, uA0, vA_Des, vA_Des_dot, FlagAttInObs, BetaAv0, SpeedA0, F_A, F_A_dot] = controlAttacker( XA, XA_goal, vA0, FlagAttInObs,BetaAv0,SpeedA0, XD, WA, NA, ND )
% ------------------------------
% Author: Vishnu S. Chipade
% Email: vishnuc@umich.edu
% ------------------------------

%Control action for a swarm of attackers, in some configuration given by
%theor communication graph and inter-agent distance potential,
%to reach protected area while avoiding obstacles and defenders
global rO rP C_d nO aO bO
global u_maxA v_maxA rho_c_A
global R_m_AO R_bar_AO R_u_AO R_m_AD R_bar_AD R_u_AD R_m_AA R_bar_AA R_u_AA
global A_A_O B_A_O C_A_O D_A_O A_A_D B_A_D C_A_D D_A_D A_A_A B_A_A C_A_A D_A_A
global E_m_AO E_bar_AO E_u_AO E_v_AO
global sigmaProdA_dot_fun
NO=size(rO,2);
global Rii00 Rij0
Rij0=R_bar_AD*ones(1,ND);
rD=XD(1:2,:);
vD=XD(3:4,:);
%Control actions for the leaders
for i=1:NA
    rA=XA(1:2,i);
    vA=XA(3:4,i);
    rA_goal=XA_goal(1:2,i);
    vA_goal=XA_goal(3:4,i);
    F_A(:,i)=[0;0];
    F_A_dot(:,i)=[0;0];
    sigmaProd=1; %for product
    sigmaBarProd=1;
    sigmaSum=0;
    sigmaBarSum=0;
    Sigma=zeros(NO+ND+NA,1);
    Sigma_dot=Sigma;
    
    % check for nearby rectangular obstacles
    if(1)
        for k=1:NO
            n=nO(k);
            a=aO(k);
            b=bO(k);
            R_OA=norm(rA-rO(1:2,k));
            R_AgoalA=norm(rA_goal-rA);
            EmAO=E_m_AO(i,k);
            EbarAO=E_bar_AO(i,k);
            EuAO=E_u_AO(i,k);
            EvAO=E_u_AO(i,k);
            E_AOk=(abs((rA(1)-rO(1,k))/a))^(2*n)+(abs((rA(2)-rO(2,k))/b))^(2*n)-1;
            E_AOk_dot=2*n*[(abs((rA(1)-rO(1,k))/a))^(2*n)/(rA(1)-rO(1,k)), (abs((rA(2)-rO(2,k))/b))^(2*n)/(rA(2)-rO(2,k))]*vA;
            %Modify the boundary elliptic distances
            E_Ades_Ok=(abs((rA_goal(1)-rO(1,k))/a))^(2*n)+(abs((rA_goal(2)-rO(2,k))/b))^(2*n)-1;
            if E_Ades_Ok<EuAO
                EuAO=E_Ades_Ok;
                EbarAO=(E_Ades_Ok-EmAO)/2;
            end
            arr_EAO(i,k)=E_AOk;
            if E_AOk>EmAO &&  E_AOk<EbarAO
                sigma=1;
                sigma_dot=0;
                sigma_bar=1;
                FlagAttInObs(i,k)=1;
            elseif  E_AOk>EbarAO &&  E_AOk<EuAO
                sigma=A_A_O(i,k)*E_AOk^3+B_A_O(i,k)*E_AOk^2+C_A_O(i,k)*E_AOk+D_A_O(i,k);
                sigma_dot=(3*A_A_O(i,k)*E_AOk^2+2*B_A_O(i,k)*E_AOk+C_A_O(i,k))*E_AOk_dot;
                sigma_bar=1;
                if FlagAttInObs(i,k)~=1
                    BetaAv0(i,k)=atan2(vA(2),vA(1)); %Orientation while entering the obstacle field
                    SpeedA0(i,k)=norm(vA); %Speed while entering the obstacle field
                end
                FlagAttInObs(i,k)=1;
            elseif E_AOk>EuAO &&  E_AOk<EvAO
                sigma_bar=A_bar_A_O(i,k)*E_AOk^3+B_bar_A_O(i,k)*E_AOk^2+C_bar_A_O(i,k)*E_AOk+D_bar_A_O(i,k);
                sigma=0;
                sigma_dot=0;
                FlagAttInObs(i,k)=0;
            else
                sigma=0;
                sigma_dot=0;
                sigma_bar=0;
                FlagAttInObs(i,k)=0;
            end
            Sigma(k)=sigma;
            Sigma_dot(k)=sigma_dot;
            
            %find vector field direction
            betaOS=atan2(rA_goal(2)-rO(2,k),rA_goal(1)-rO(1,k));  %angle between the desired location and the obstacle
            betaOA=atan2(rA(2)-rO(2,k),rA(1)-rO(1,k));   %angle between the location and the obstacle
            betaOA_dot=[-(rA(2)-rO(2,k)), (rA(1)-rO(1,k))]*vA/R_OA^2;
            beta_barOS=atan2(b^(2*n)*sign(rA_goal(1)-rO(1,k))*(abs(rA_goal(1)-rO(1,k)))^(2*n-1),-a^(2*n)*sign(rA_goal(2)-rO(2,k))*(abs(rA_goal(2)-rO(2,k)))^(2*n-1));  %angle of the tangent at the safe area
            betaSA=atan2(rA(2)-rA_goal(2),rA(1)-rA_goal(1));
            betaSA_dot=[-(rA(2)-rA_goal(2)), (rA(1)-rA_goal(1))]*(vA-vA_goal)/R_AgoalA^2;
            beta_bar=atan2(b^(2*n)*sign(rA(1)-rO(1,k))*(abs(rA(1)-rO(1,k)))^(2*n-1),-a^(2*n)*sign(rA(2)-rO(2,k))*(abs(rA(2)-rO(2,k)))^(2*n-1));  %angle of the tangent at the given location
            if beta_bar<0   % get betaOS between [0,2pi]
                beta_bar=beta_bar+2*pi;
            end
            Delta_beta=betaOA-BetaAv0(i,k);
            if Delta_beta<0
                Delta_beta=Delta_beta+2*pi;
            end
            
            Delta_beta_bar=beta_barOS-betaOS;
            if Delta_beta_bar<0
                Delta_beta_bar=Delta_beta_bar+2*pi;
            end
            
            crossProd=cross([cos(BetaAv0(i,k));sin(BetaAv0(i,k));0],[cos(beta_bar);sin(beta_bar);0]);
            if crossProd(3)>0
                beta_bar1=BetaAv0(i,k);
                beta_bar1_dot=0;
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
                beta_bar1_dot=(b^(2*n)*(2*n-1)*sign(cos(betaOA))*abs(cos(betaOA))^(2*n-2))/(a^(2*n)*sign(sin(betaOA))*abs(sin(betaOA))^(2*n))*betaOA_dot;
            end
            F_AOk=  [cos(beta_bar1);sin(beta_bar1)];
            F_A(:,i) = F_A(:,i) + sigma*F_AOk;
            F_A_dot(:,i)=F_A_dot(:,i)+sigma_dot*F_AOk+sigma*beta_bar1_dot*[-sin(beta_bar1);cos(beta_bar1)];
            sigmaProd=sigmaProd*(1-sigma);
            sigmaSum=sigmaSum+sigma;
            sigmaBarProd=sigmaBarProd*(1-sigma_bar);
            sigmaBarSum=sigmaBarSum+sigma_bar;
        end
    end
    
    sigmaSumA=0;
    sigmaProdA=1;
    %check for nearby attackers
    if(1)
        for ii=1:NA
            if ii~=i
                rAAii=XA(1:2,ii)-rA;
                R_AAii=norm(rAAii);
                R_AAii_dot=(-rAAii)'*vA/R_AAii;
                if R_AAii<=rho_c_A
                    if R_AAii>R_m_AA && R_AAii<R_bar_AA
                        sigma=1;     sigma_dot=0;
                    elseif R_AAii>R_bar_AA && R_AAii<R_u_AA
                        sigma=A_A_A*R_AAii^3+B_A_A*R_AAii^2+C_A_A*R_AAii+D_A_A;
                        sigma_dot=(3*A_A_A*R_AAii^2+2*B_A_A*R_AAii+C_A_A)*R_AAii_dot;
                    else
                        sigma=0;   sigma_dot=0;
                    end
                    Sigma(NO+ND+ii)=sigma;
                    Sigma_dot(NO+ND+ii)=sigma_dot;
                    F_AAii=(-rAAii)/R_AAii;
                    F_A(:,i) = F_A(:,i) + sigma*F_AAii;
                    F_A_dot(:,i) = F_A_dot(:,i) + sigma_dot*F_AAii+sigma*(vA-R_AAii_dot*F_AAii)/R_AAii;
                    sigmaProdA=sigmaProdA*(1-sigma);
                    sigmaSumA=sigmaSumA+sigma;
                end
            end
        end
    end
    
    %for attraction toward the protected area
    rAP=XA_goal(1:2,i)-rA;
    R_AP=norm(rAP);
    R_AP_dot=(rAP)'*(-vA)/R_AP;
    F_AP=(rAP)/R_AP;
    F_A(:,i) = F_A(:,i) + sigmaProd*F_AP;%/Rc^2;
    s1=num2cell([Sigma',Sigma_dot']);
    sigmaProd_dot=sigmaProdA_dot_fun(s1{:});
    F_A_dot(:,i)=F_A_dot(:,i)+sigmaProd_dot*F_AP+sigmaProd*(-vA-R_AP_dot*F_AP)/R_AP;
    norm_F_Ai=norm(F_A(:,i));
    if norm_F_Ai~=0
        vA_hat=F_A(:,i)/norm_F_Ai;
    else
        vA_hat=[-vA(2,i),vA(1)]'/norm(vA(1:2,i));
    end
    % vA_des=.75*v_maxA(i)*tanh(R_AP)*vA_hat;
    indAttInObs=find(FlagAttInObs(i,:)==1);
    if ~isempty(indAttInObs)
        vA_des=SpeedA0(i,indAttInObs)*vA_hat;
        norm_F_Ai_dot=F_A(:,i)'*F_A_dot(:,i)/norm_F_Ai;
        vA_hat_dot=(F_A_dot(:,i)-vA_hat*norm_F_Ai_dot)/norm_F_Ai;
        %    vA_des_dot=0.75*v_maxA(i)*tanh(R_AP)*vA_hat_dot+0.75*v_maxA(i)*(1-tanh(R_AP)^2)*R_AP_dot*vA_hat;
        vA_des_dot=SpeedA0(i,indAttInObs)*vA_hat_dot;
    else
        vA_des=zeros(2,1);
        vA_des_dot=vA_des;
    end
    %     uA(:,i)=-(vA-vA_des)+vA_des_dot+C_d*vA_des;
    vA_Des(:,i)=vA_des;
    vA_Des_dot(:,i)=vA_des_dot;
    
    
    %check formation controller
    %Vij(x)=log(k/(x-x0)+(x-x0)/k)
    nabla_ri_Vii1=zeros(2,1);
    vNA=zeros(2,1);
    for ii=find(WA(i,:)==1)
        Rii0=Rii00(i);
        Ri_ii=norm(rA-XA(1:2,ii));
        nabla_ri_Vii1=nabla_ri_Vii1+WA(i,ii)*(rA-XA(1:2,ii))/Ri_ii^2*((Ri_ii-R_bar_AA)^2-Rii0^2)/((Ri_ii-R_bar_AA)^2+Rii0^2);
        %nabla_ri_Vii1=nabla_ri_Vii1+WA(i,ii)*(rA-XA(1:2,ii))*((Ri_ii)^2-Rii0^2);
        vNA=vNA+WA(i,ii)*XA(3:4,i);
    end
    
    
    %For defenders nearby
    nabla_ri_Vij=zeros(2,1);
    for j=1:ND
        Rij=Rij0(j);
        Ri_j=norm(rA-rD(1:2,j));
        nabla_ri_Vij=nabla_ri_Vij-(rA-rD(1:2,j))/Ri_j^2*((Ri_j-R_bar_AD)^2-Rij^2)/((Ri_j-R_bar_AD)^2+Rij^2);
        %vNA=vNA+WA(i,ii)*XA(3:4,i);
    end
    
    
    kf=5;
    if sum(WA(i,:))~=0
        uA00(:,i)=-0.9*(vA-vNA/sum(WA(i,:)))-kf*nabla_ri_Vii1;%0.01*(vA-vA_des)+C_d*vA_des;  %add vA_des_dot
    else
        uA00(:,i)=-kf*nabla_ri_Vii1;%;-(vA-vA_des)+C_d*vA_des;
    end
    
    %check for nearby defenders
    F_AD(:,i)=[0,0]';
    F_AD_dot(:,i)=[0,0]';
    sigmaSumD=0;
    sigmaProdD=1;
    minRAD=Inf;
    uAD=zeros(2,1);
    if(1)
        for j=1:ND
            rADj=rD(1:2,j)-rA;
            R_ADj=norm(rADj);
            R_ADj_dot=(-rADj)'*vA/R_ADj;
            if R_ADj<=rho_c_A
                if R_ADj<R_u_AD
                    if R_ADj>R_m_AD && R_ADj<R_bar_AD
                        sigma=1;     sigma_dot=0;
                    elseif R_ADj>R_bar_AD && R_ADj<R_u_AD
                        sigma=A_A_D*R_ADj^3+B_A_D*R_ADj^2+C_A_D*R_ADj+D_A_D;
                        sigma_dot=(3*A_A_D*R_ADj^2+2*B_A_D*R_ADj+C_A_D)*R_ADj_dot;
                    end
                    nabla_ri_Vij=-(rA-rD(1:2,j))/Ri_j^2*((Ri_j-R_bar_AD)^2-Rij^2)/((Ri_j-R_bar_AD)^2+Rij^2);
                    uAD=uAD-sigma*(vA-vD(:,j))-sigma*nabla_ri_Vij;
                else
                    sigma=0;   sigma_dot=0;
                end
                Sigma(NO+j)=sigma;
                Sigma_dot(NO+j)=sigma_dot;
                F_ADj=(-rADj)/R_ADj;
                F_AD(:,i) = F_AD(:,i) + sigma*F_ADj;
                F_AD_dot(:,i) = F_AD_dot(:,i) + sigma_dot*F_ADj+sigma*(vA-R_ADj_dot*F_ADj)/R_ADj;
                sigmaProdD=sigmaProdD*(1-sigma);
                sigmaSumD=sigmaSumD+sigma;
            end
            if minRAD>R_ADj
                minRAD=R_ADj;
            end
        end
    end
    
    
    duA(:,i)=  1*F_AD(:,i)*(R_u_AD-minRAD)/(minRAD-R_m_AD);
    %duA(:,i)=  F_AD(:,i)*0.5*(tanh(-4/(R_u_AD-R_m_AD)*(minRAD-R_m_AD)-2)+1);
    uA0(:,i)=sigmaProd*sigmaProdD*(uA00(:,i)-0.01*(rA-rA_goal));
    uA(:,i)=uA0(:,i)+sigmaSumD*duA(:,i)+sigmaSum*(-(vA-vA_des)-sigmaSumA*kf*nabla_ri_Vii1);
    uA(:,i)=uA0(:,i)+uAD+sigmaSum*(-(vA-vA_des)-sigmaSumA*kf*nabla_ri_Vii1);
    
    infNorm_uA=max(abs(uA(1,i)),abs(uA(2,i)));
    if infNorm_uA>u_maxA(i)
        uA(:,i)=uA(:,i)*u_maxA(i)/infNorm_uA;
    end
    
    % uA(:,i)=uA00(:,i);
    vA_dot(:,i)=uA(:,i);
    
end

%Control actions for the followers


end
