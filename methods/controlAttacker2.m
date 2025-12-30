function [ uA, uA0, vA_Des, vA_Des_dot, rAProj1, FlagAttInObs, BetaAv0, SpeedA0, SigmaProdD,F_A, F_A_dot] = controlAttacker2( XA, XA_goal, vA0, FlagAttInObs,BetaAv0,SpeedA0, XD, WA, NA, ND )
% ------------------------------
% Author: Vishnu S. Chipade
% Email: vishnuc@umich.edu
% ------------------------------

%Control action for a swarm of attackers, in some configuration given by
%theor communication graph and inter-agent distance potential,
%to reach protected area while avoiding obstacles and defenders

%Obstacle avoidance using beta-agent idea rather than vector fields for the
%flock
global rO rP C_d nO aO bO
global u_maxA v_maxA rho_c_A
global R_m_AO R_bar_AO R_u_AO R_m_AD R_bar_AD R_u_AD R_m_AA R_bar_AA R_u_AA
global A_A_O B_A_O C_A_O D_A_O A_A_D B_A_D C_A_D D_A_D A_A_A B_A_A C_A_A D_A_A
global E_m_AO E_bar_AO E_u_AO E_v_AO
global sigmaProdA_dot_fun
global A_bar_A_O B_bar_A_O C_bar_A_O D_bar_A_O
NO=size(rO,2);
global Rii00 Rij0 Rik00 options
global kd

rD=XD(1:2,:);
vD=XD(3:4,:);


rAcm=sum(XA(1:2,:),2)/NA;
vAcm=sum(XA(3:4,:),2)/NA;

%Control actions for the leaders
SigmaProdD=ones(NA,1);
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
            rAProj(:,k)=rO(1:2,k);
            vAProj(:,k)=[0,0]';
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
            if E_AOk<EuAO
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
                end
                %Point along the radial direction as the initial condition 
                beta=atan2(rA(2)-rO(2,k),rA(1)-rO(1,k));
                RProj=((EmAO+1)/(abs(cos(beta))^(2*n)/a^(2*n)+abs(sin(beta))^(2*n)/b^(2*n)))^(1/(2*n));
                rP0=[RProj*cos(beta);RProj*sin(beta)]+rO(1:2,k);
                f=@(r) [abs(r(1)-rO(1,k))^(2*n)/a^(2*n)+abs(r(2)-rO(2,k))^(2*n)/b^(2*n)-EmAO-1; a^(2*n)*sign(r(2)-rO(2,k))*abs(r(2)-rO(2,k))^(2*n-1)*(rA(1)-r(1))-b^(2*n)*sign(r(1)-rO(1,k))*abs(r(1)-rO(1,k))^(2*n-1)*(rA(2)-r(2))];
                rAProj(:,k)=fsolve(f,rP0,options);
                %Tangent at the projection point
                beta_barP=atan2(b^(2*n)*sign(rAProj(1,k)-rO(1,k))*abs(rAProj(1,k)-rO(1,k))^(2*n-1),-a^(2*n)*sign(rAProj(2,k)-rO(2,k))*abs(rAProj(2,k)-rO(2,k))^(2*n-1));
                rTP=[cos(beta_barP);sin(beta_barP)];
                if rTP'*vA<0  %If the velocity of A is opposite to the tangent then reveresed the tangent direction for projection
                    rTP=-rTP;
                end
                vAProj(:,k)=rTP'*vA*rTP;
                
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
            
            sigmaProd=sigmaProd*(1-sigma);
            sigmaSum=sigmaSum+sigma;
            sigmaBarProd=sigmaBarProd*(1-sigma_bar);
            sigmaBarSum=sigmaBarSum+sigma_bar;
        end
    end
       
    
    %check formation controller
    %Vij(x)=log(k/(x-x0)+(x-x0)/k)
    nabla_ri_Vii=zeros(2,1);
    vNA=zeros(2,1);
    for ii=find(WA(i,:)==1)
        Rii0=Rii00(i);
        Ri_ii=norm(rA-XA(1:2,ii));
        nabla_ri_Vii=nabla_ri_Vii+WA(i,ii)*(rA-XA(1:2,ii))/Ri_ii*(Ri_ii-R_m_AA)*((Ri_ii-R_m_AA)^2-Rii0^2)/((Ri_ii-R_m_AA)^2+Rii0^2);
        %nabla_ri_Vii=nabla_ri_Vii+WA(i,ii)*(rA-XA(1:2,ii))*((Ri_ii)^2-Rii0^2);
        vNA=vNA+WA(i,ii)*XA(3:4,i);
    end
    
    %Potential for projected beta-agents on the static obstacles
    ko=.1;
    nabla_ri_ViP=zeros(2,1);
    vNP=zeros(2,1);
    for k=1:NO
        if Sigma(k)~=0
        Rik0=Rik00(i);
        Ri_ik=norm(rA-rAProj(:,k));
        nabla_ri_ViP= nabla_ri_ViP+Sigma(k)*(rA-rAProj(:,k))/Ri_ik/(Ri_ik-R_m_AA)*((Ri_ik-R_m_AA)^2-Rik0^2)/((Ri_ik-R_m_AA)^2+Rik0^2);
        vNP=vNP+Sigma(k)*vAProj(:,k);
        end
    end
    if sum(Sigma)~=0
        uAP0(:,i)=-1*(sum(Sigma)*vA-vNP)-ko*nabla_ri_ViP;%0.01*(vA-vA_des)+C_d*vA_des;  %add vA_des_dot
    else
        uAP0(:,i)=-ko*nabla_ri_ViP;%;-(vA-vA_des)+C_d*vA_des;
    end
    
%      %For defenders nearby
%     nabla_ri_Vij=zeros(2,1);
%     for j=1:ND
%         Rij=Rij0(j);
%         Ri_j=norm(rA-rD(1:2,j));
%         nabla_ri_Vij=nabla_ri_Vij-(rA-rD(1:2,j))/Ri_j^2*((Ri_j-R_bar_AD)^2-Rij^2)/((Ri_j-R_bar_AD)^2+Rij^2);
%         %vNA=vNA+WA(i,ii)*XA(3:4,i);
%     end
    
    kf=.1;
    if sum(WA(i,:))~=0
        uA00(:,i)=-1*(sum(WA(i,:))*vA-vNA)-kf*nabla_ri_Vii;%0.01*(vA-vA_des)+C_d*vA_des;  %add vA_des_dot
    else
        uA00(:,i)=-kf*nabla_ri_Vii;%;-(vA-vA_des)+C_d*vA_des;
    end
    
    %check for nearby defenders
    F_AD(:,i)=[0,0]';
    F_AD_dot(:,i)=[0,0]';
    sigmaSumD=0;
    sigmaProdD=1;
    minRAD=Inf;
    uADv=zeros(2,1);
    uADr=zeros(2,1);
    if(1)
        for j=1:ND
            rADj=rD(1:2,j)-rA;
            R_ADj=norm(rADj);
            R_ADj_dot=(-rADj)'*vA/R_ADj;
            if R_ADj<=rho_c_A
                if R_ADj<R_u_AD
                    if R_ADj>R_m_AD && R_ADj<R_bar_AD
                        sigma=1;     sigma_dot=0;
                        ind=find(WA(i,:)==1);
                        SigmaProdD(ind)=SigmaProdD(ind)*(1-sigma);
                    elseif R_ADj>R_bar_AD && R_ADj<R_u_AD
                        sigma=A_A_D*R_ADj^3+B_A_D*R_ADj^2+C_A_D*R_ADj+D_A_D;
                        sigma_dot=(3*A_A_D*R_ADj^2+2*B_A_D*R_ADj+C_A_D)*R_ADj_dot;
                    end
                    Rij=Rij0(j);
                    nabla_ri_Vij=kd*(rA-rD(1:2,j))/R_ADj/(R_ADj-R_m_AD)*((R_ADj-R_m_AD)^2-Rij^2)/((R_ADj-R_m_AD)^2+Rij^2);
                    uADv=uADv-sigma*(vA-vD(:,j));
                    uADr=uADr-sigma*nabla_ri_Vij;
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
    
    %uA(:,i)=uA0(:,i)+sigmaSumD*duA(:,i)+uAP0(:,i);
    
    if (1)
       %uA0(:,i)= uA00(:,i)-sigmaProdD*(0.005*(rA-rA_goal)+0.5*(vA-vA_goal))+uAP0(:,i);
    uA(:,i)=uA00(:,i)+uADr-sigmaProdD*(0.005*(rA-rA_goal)+0.5*(vA-vA_goal))+uAP0(:,i); %+uADv
    uA0(:,i)=uA(:,i);
   % uA(:,i)=-sigmaProdD*(0.005*(rA-rA_goal)+0.5*(vA-vA_goal));
    else
    uA(:,i)=uA00(:,i)+uADv+uADr+uAP0(:,i);
    end
    
    infNorm_uA=max(abs(uA(1,i)),abs(uA(2,i)));
    if infNorm_uA>u_maxA(i)
        uA(:,i)=uA(:,i)*u_maxA(i)/infNorm_uA;
    end
    
    % uA(:,i)=uA00(:,i);
    vA_dot(:,i)=uA(:,i);
    vA_Des(:,i)=zeros(2,1);
    vA_Des_dot(:,i)=zeros(2,1);
    rAProj1(2*i-1:2*i,:)=rAProj;
    SigmaProdD(i,1)=sigmaProdD;
end

%Control actions for the followers


end
