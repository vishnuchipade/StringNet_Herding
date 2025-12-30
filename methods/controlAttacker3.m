function [ uA, uA0, R_AO_min, R_AAProjS_min, vA_Des, vA_Des_dot, rAProj1, FlagAttInObs, BetaAv0, SpeedA0,mAv0,cAv0, SigmaProdD,F_A, F_A_dot] = controlAttacker3( XA, XA_goal, vA0, FlagAttInObs,BetaAv0,SpeedA0,mAv0,cAv0, XD, WA, WDString, NA, ND )
% ------------------------------
% Author: Vishnu S. Chipade
% Email: vishnuc@umich.edu
% ------------------------------

%Control action for a swarm of attackers, in some configuration given by
%their communication graph and inter-agent distance potential,
%to reach protected area while avoiding obstacles and defenders

%Obstacle avoidance using beta-agent idea rather than vector fields for the
%flock
global obs rP C_d rhoAD_safe
global u_maxA v_maxA rho_c_A
global R_m_AO R_bar_AO R_u_AO R_m_AD R_bar_AD R_u_AD R_m_AA R_bar_AA R_u_AA
global A_A_O B_A_O C_A_O D_A_O A_A_D B_A_D C_A_D D_A_D A_A_A B_A_A C_A_A D_A_A
global sigmaProdA_dot_fun
global Rii_bar Rij0 Rik00 options
global kd largeP rho_Acon
global kADr kAFr kAOr kAOr2 kAPr kAPv kAOv2 kAFv alphaAFv kAOv alphaAOv kADv alphaADv

NO=obs.NO;
rCO2=obs.rCO2;
rVO=obs.rVO1;  %to have the first vertex repeated
if ND>0
    rD=XD(1:2,:);
    vD=XD(3:4,:);
end

rAcm=sum(XA(1:2,:),2)/NA;
vAcm=sum(XA(3:4,:),2)/NA;

%Control actions
SigmaProdD=ones(NA,1);
R_AO_min=Inf;
R_AAProjS_min=Inf;  %distance from the strings
for i=1:NA
    WDString_temp=WDString;
    rA=XA(1:2,i);
    vA=XA(3:4,i);
    if norm(vA)>v_maxA(i)
        vA=vA*v_maxA(i)/norm(vA);
    end
    
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
    
    % check for collision avoidance from nearby obstacles
    %Potential for projected beta-agents on the static obstacles
    uAOv=zeros(2,1);
    uAOr=uAOv;
    if(1)
        for k=1:NO
            R_AOk=Inf;
            %rVO{k}=[rVO{k},rVO{k}(:,1)];
            for kk=1:length(rVO{k}(1,:))-1
                dx=rVO{k}(1,kk+1)-rVO{k}(1,kk);
                dy=rVO{k}(2,kk+1)-rVO{k}(2,kk);
                if abs(dx)>1e-10
                    mVV=dy/dx;  %slope of the line
                    cVV=rVO{k}(2,kk)-mVV*rVO{k}(1,kk);
                    rAProjO_temp(1,1)=(mVV*rA(2)+rA(1)-mVV*cVV)/(1+mVV^2);
                    rAProjO_temp(2,1)=mVV*rAProjO_temp(1,1)+cVV;
                    lambdaDP=(rAProjO_temp(1,1)-rVO{k}(1,kk))/dx;
                else
                    rAProjO_temp(1,1)=rVO{k}(1,kk);
                    rAProjO_temp(2,1)=rA(2);
                    lambdaDP=(rAProjO_temp(2,1)-rVO{k}(2,kk))/dy;
                end
                if lambdaDP<0
                    rAProjO_temp(:,1)=rVO{k}(1:2,kk);
                elseif lambdaDP>1
                    rAProjO_temp(:,1)=rVO{k}(1:2,kk+1);
                end
                R_AOk_temp=norm(rA-rAProjO_temp);
                if R_AOk_temp<R_AOk
                    R_AOk=R_AOk_temp;
                    rAProjO(:,k)=rAProjO_temp;
                end
            end
            
            R_AOk=R_AOk-rhoAD_safe;
            if R_AOk<R_AO_min
                R_AO_min=R_AOk;
            end
            theta=atan2(rA(2)-rAProjO(2,k),rA(1)-rAProjO(1,k));
            rTP=[cos(theta+pi/2), sin(theta+pi/2)]';
            if rTP'*vA<0  %If the velocity of D is opposite to the tangent then reveres the tangent direction for projection
                rTP=-rTP;
            end
            vAProjO(:,k)=rTP'*vA*rTP;
            %arr_RAO(i,k)=R_AOk;
            if R_AOk<R_u_AO(k)
                if R_AOk<R_bar_AO(k)  %&& E_AOk>EmAO
                    sigma=1;
                    sigma_dot=0;
                    sigma_bar=1;
                    FlagAttInObs(i,k)=1;
                elseif  R_AOk>R_bar_AO(k) &&  R_AOk<R_u_AO(k)
                    sigma=A_A_O(i,k)*R_AOk^3+B_A_O(i,k)*R_AOk^2+C_A_O(i,k)*R_AOk+D_A_O(i,k);
                    %sigma_dot=(3*A_A_O(i,k)*R_AOk^2+2*B_A_O(i,k)*R_AOk+C_A_O(i,k))*R_AOk_dot;
                    sigma_bar=1;
                end
            else
                sigma=0;
                sigma_dot=0;
                sigma_bar=0;
                FlagAttInObs(i,k)=0;
            end
            Sigma(k)=sigma;
            %Sigma_dot(k)=sigma_dot;
            
            sigmaProd=sigmaProd*(1-sigma);
            sigmaSum=sigmaSum+sigma;
            sigmaBarProd=sigmaBarProd*(1-sigma_bar);
            sigmaBarSum=sigmaBarSum+sigma_bar;
            
            if Sigma(k)~=0
                Rik0=Rik00(i);
                %find the projection of the obstacle center on the line joining
                %current position and goal position to check if the obstacle is
                %in between or not
                dx=rA(1)-rA_goal(1);
                dy=rA(2)-rA_goal(2);
                mAAG=dy/dx;  %slope of line joining the defenders
                cAAG=rA(2)-mAAG*rA(1);
                rOProj(1)=(mAAG*rCO2(2,k)+rCO2(1,k)-mAAG*cAAG)/(1+mAAG^2);
                rOProj(2)=mAAG*rOProj(1)+cAAG;
                if mAAG<1e16
                    lambdaOP=(rOProj(1)-rA_goal(1))/dx;
                else
                    lambdaOP=(rOProj(2)-rA_goal(2))/dy;
                end
                
                if (R_AOk-R_m_AA)>1e-1
                    nabla_ri_ViP= kAOr*Sigma(k)*(rA-rAProjO(:,k))/R_AOk/abs(R_AOk-R_m_AA)*((R_AOk-R_m_AA)^2-Rik0^2)/((R_AOk-R_m_AA)^2+Rik0^2);
                else
                    nabla_ri_ViP= - kAOr*Sigma(k)*(rA-rAProjO(:,k))*largeP;
                end
                if norm(vA-vAProjO(:,k))>1e-16
                    dv=kAOv*(vA-vAProjO(:,k))*norm(vA-vAProjO(:,k))^(alphaAOv-1);
                else
                    dv=[0,0]';
                end
                
                uAOr=uAOr-Sigma(k)*nabla_ri_ViP;
                if lambdaOP<=1 && lambdaOP>=0  %if the obstacle is in the way
                    uAOv=uAOv-Sigma(k)*dv;
                else  %Obstacle is no longer in way
                    uAOv=uAOv+[0,0]';
                end
            end
        end
    end
    
    
    %check formation controller
    %Vij(x)=log(k/(x-x0)+(x-x0)/k)
    
    uAFv=zeros(2,1);
    uAFr=uAFv;
    for ii=find(WA(i,:)==1)
        Rii0=Rii_bar(i,ii)-R_m_AA;  %k
        Ri_ii=norm(rA-XA(1:2,ii));
        if abs(Ri_ii-R_m_AA)>1e-1
            nabla_ri_Vii=kAFr*WA(i,ii)*(rA-XA(1:2,ii))/Ri_ii/abs(Ri_ii-R_m_AA)*((Ri_ii-R_m_AA)^2-Rii0^2)/((Ri_ii-R_m_AA)^2+Rii0^2);
        else
            nabla_ri_Vii=-kAFr*WA(i,ii)*(rA-XA(1:2,ii))*largeP;
        end
        if norm(vA-XA(3:4,ii))>1e-16
            dv=kAFv*(vA-XA(3:4,ii))*norm(vA-XA(3:4,ii))^(alphaAFv-1);
        else
            dv=[0,0]';
        end
        uAFv=uAFv-dv;
        uAFr=uAFr-nabla_ri_Vii;
    end
    
    %check for nearby defenders
    F_AD(:,i)=[0,0]';
    F_AD_dot(:,i)=[0,0]';
    sigmaSumD=0;
    sigmaProdD=1;
    minRAD=Inf;
    uADv=zeros(2,1);
    uADr=zeros(2,1);
    countAPS=0;
    if(1)
        if i>15
            [uAD_pot, minRAD]=potentialControl(XA(:,i),XD,5*R_m_AD,2*R_bar_AD, 2*R_u_AD, 2*Rij0(1),kADr,kADv, alphaADv);
        else
             [uAD_pot, minRAD]=potentialControl(XA(:,i),XD,R_m_AD,R_bar_AD, R_u_AD, Rij0(1),kADr,kADv, alphaADv);
        end
        for j=1:ND
%             rADj=rD(1:2,j)-rA;
%             R_ADj=norm(rADj);
%             R_ADj_dot=(-rADj)'*vA/R_ADj;
%             if R_ADj<=rho_c_A
%                 if R_ADj<R_u_AD
%                     if R_ADj>R_m_AD && R_ADj<R_bar_AD
%                         sigma=1;     sigma_dot=0;
%                         ind=find(WA(i,:)==1);
%                         SigmaProdD(ind)=SigmaProdD(ind)*(1-sigma);
%                     elseif R_ADj>R_bar_AD && R_ADj<R_u_AD
%                         sigma=A_A_D*R_ADj^3+B_A_D*R_ADj^2+C_A_D*R_ADj+D_A_D;
%                         sigma_dot=(3*A_A_D*R_ADj^2+2*B_A_D*R_ADj+C_A_D)*R_ADj_dot;
%                     end
%                     Rij=Rij0(j);
%                     if (R_ADj-R_m_AD)>1e-1
%                         nabla_ri_Vij=kADr*(rA-rD(1:2,j))/R_ADj/abs(R_ADj-R_m_AD)*((R_ADj-R_m_AD)^2-Rij^2)/((R_ADj-R_m_AD)^2+Rij^2);
%                     else
%                         nabla_ri_Vij=-kADr*(rA-rD(1:2,j))*largeP;
%                     end
%                     dv=kADv*(vA-vD(:,j))*norm(vA-vD(:,j))^(alphaADv-1);
%                     uADv=uADv-sigma*dv;
%                     uADr=uADr-sigma*nabla_ri_Vij;
%                 else
%                     sigma=0;   sigma_dot=0;
%                 end
%                 Sigma(NO+j)=sigma;
%                 Sigma_dot(NO+j)=sigma_dot;
%                 F_ADj=(-rADj)/R_ADj;
%                 F_AD(:,i) = F_AD(:,i) + sigma*F_ADj;
%                 F_AD_dot(:,i) = F_AD_dot(:,i) + sigma_dot*F_ADj+sigma*(vA-R_ADj_dot*F_ADj)/R_ADj;
%                 sigmaProdD=sigmaProdD*(1-sigma);
%                 sigmaSumD=sigmaSumD+sigma;
%             end
%             if minRAD>R_ADj
%                 minRAD=R_ADj;
            % end
            
            
            %Find the projections of the attacker on the string
            %attached to defender j
            for jj=find(WDString_temp(j,:)==1)
                dx=XD(1,jj)-XD(1,j);
                dy=XD(2,jj)-XD(2,j);
                mVV=dy/dx;  %slope of line joining the defenders
                cVV=XD(2,j)-mVV*XD(1,j);
                countAPS=countAPS+1;
                rAProjS(1,countAPS)=(mVV*rA(2)+rA(1)-mVV*cVV)/(1+mVV^2);
                rAProjS(2,countAPS)=mVV*rAProjS(1,countAPS)+cVV;
                if mVV<1e16
                    lambdaAP=(rAProjS(1,countAPS)-XD(1,j))/dx;
                else
                    lambdaAP=(rAProjS(2,countAPS)-XD(2,j))/dy;
                end
                if lambdaAP<0
                    rAProjS(:,countAPS)=XD(1:2,j);
                elseif lambdaAP>1
                    rAProjS(:,countAPS)=XD(1:2,jj);
                end
                rTDD=XD(1:2,j)-XD(1:2,jj);
                rTDD=rTDD/norm(rTDD);
                if rTDD'*vA<0
                    rTDD=-rTDD;
                end
                %velocity projection due to Attacker
                vAProjS0(:,countAPS)=rTDD'*vA*rTDD;
                %velocity because of the motion of the defenders
                vAProjS1(:,countAPS)=XD(3:4,j)+(XD(3:4,jj)-XD(3:4,j))*norm(rAProjS(1:2,countAPS)-XD(1:2,j))/norm(XD(1:2,jj)-XD(1:2,j));
                vAProjS(:,countAPS)=vAProjS0(:,countAPS)+vAProjS1(:,countAPS);
                %Remove the strings from calculations
                WDString_temp(j,jj)=0;
                WDString_temp(jj,j)=0;
            end
        end
    end
    
    %Check for the nearby strings
    uAAProjS=zeros(2,1);
    uAAProjSv=zeros(2,1);
    sigmaProdDS=1;
    for js=1:countAPS
        rAAProjS=rAProjS(1:2,js)-rA;
        R_AAProjS=norm(rAAProjS);
        if R_AAProjS<R_AAProjS_min
            R_AAProjS_min=R_AAProjS;
        end
        if R_AAProjS<R_u_AD
            if R_AAProjS<R_bar_AD
                sigma=1;     sigma_dot=0;
            elseif R_AAProjS>R_bar_AD && R_AAProjS<R_u_AD
                sigma=A_A_D*R_AAProjS^3+B_A_D*R_AAProjS^2+C_A_D*R_AAProjS+D_A_D;
                %sigma_dot=(3*A_A_D*R_AAProjS^2+2*B_A_D*R_AAProjS+C_A_D)*R_AAProjS_dot;
            end
            Rij=Rij0(1);   %Any large value to avoid
            if (R_AAProjS-R_m_AD)>1e-1
                nabla_ri_ViPS=kADr*(-rAAProjS)/R_AAProjS/abs(R_AAProjS-R_m_AD)*((R_AAProjS-R_m_AD)^2-Rij^2)/((R_AAProjS-R_m_AD)^2+Rij^2);
            else
                nabla_ri_ViPS=kADr*(rAAProjS)*largeP;
            end
            % uADv=uADv-sigma*(vA-vD(:,j));
            uAAProjS=uAAProjS-sigma*nabla_ri_ViPS;
            if norm(vA-vAProjS(:,js))>1e-16 && ((vA-vAProjS(:,js))'*(rA-rAProjS(:,js))<0)
                dv=kADv*(vA-vAProjS(:,js))*norm(vA-vAProjS(:,js))^(alphaADv-1);
            else
                dv=[0,0]';
            end
            uAAProjSv=uAAProjSv - sigma*dv;
        else
            sigma=0;   sigma_dot=0;
        end
        sigmaProdDS=sigmaProdDS*(1-sigma);
    end
    
    
    duA(:,i)=  1*F_AD(:,i)*(R_u_AD-minRAD)/(minRAD-R_m_AD);
    %duA(:,i)=  F_AD(:,i)*0.5*(tanh(-4/(R_u_AD-R_m_AD)*(minRAD-R_m_AD)-2)+1);
    
    %uA(:,i)=uA0(:,i)+sigmaSumD*duA(:,i)+uAP0(:,i);
    
    %Stay within the connectivity region of the attackers
    if (0)
        if norm(rA-rAcm)<rho_Acon
            thetaAAcm=atan2(rA(2)-rAcm(2),rA(1)-rAcm(1));
            rAProjC=rAcm+rho_Acon*[cos(thetaAAcm),sin(thetaAAcm)]';
            rTP=[cos(thetaAAcm+pi/2),sin(thetaAAcm+pi/2)]';
            if rTP'*vA<0  %If the velocity of A is opposite to the tangent then reveres the tangent direction for projection
                rTP=-rTP;
            end
            vAProjC=rTP'*vA*rTP+vAcm;
            Rik0=Rik00(1);
            Ri_ik=norm(rA-rAProjC);
            nabla_ri_ViP= kAOr2*(rA-rAProjC)/Ri_ik/abs(Ri_ik-R_m_AA)*((Ri_ik-R_m_AA)^2-Rik0^2)/((Ri_ik-R_m_AA)^2+Rik0^2);
            if norm(vA-vAProjC)>1e-16
                dvP=-kAOv2*(vA-vAProjC)*norm(vA-vAProjC)^(alphaAOv-1);
            else
                dvP=[0,0]';
            end
        else
            nabla_ri_ViP=[0,0]';
            dvP=[0,0]';
        end
    end
    if (1)
        %sigmaProdD=0;
        %uA0(:,i)= uA00(:,i)-sigmaProdD*(0.005*(rA-rA_goal)+0.5*(vA-vA_goal))+uAP0(:,i);
        duA_goal=-(kAPr*(rA-rA_goal)+kAPv*(vA));  %sigmaProdD*sigmaProdDS*
        norm_duA_goal=norm(duA_goal);
        if norm_duA_goal>1e-10
            uA(:,i)=min(1*u_maxA(i),norm_duA_goal)*(duA_goal/norm_duA_goal)+uAFv+uAFr+uAOv+uAOr+uAD_pot+uAAProjS+uAAProjSv+C_d*norm(vA)*vA;% +dvP -nabla_ri_ViP; %+uADv
        else
            uA(:,i)=uAFv+uAFr+uAOv+uAOr+uAD_pot+uAAProjS+uAAProjSv+C_d*norm(vA)*vA;% +dvP -nabla_ri_ViP; %+uADv
        end
        uA0(:,i)=uA(:,i);
        % uA(:,i)=-sigmaProdD*(0.005*(rA-rA_goal)+0.5*(vA-vA_goal));
    else
        uA(:,i)=uA00(:,i)+uADv+uADr+uAP0(:,i);
    end
    
    %infNorm_uA=max(abs(uA(1,i)),abs(uA(2,i)));
    norm_uA=norm(uA(:,i));
    if norm_uA>u_maxA(i) %infNorm_uA>u_maxA(i)
        %uA(:,i)=uA(:,i)*u_maxA(i)/infNorm_uA;
        uA(:,i)=uA(:,i)*u_maxA(i)/norm_uA;
    end
    if isnan(uA(1,i))
        
    end
    
    % uA(:,i)=uA00(:,i);
    vA_dot(:,i)=uA(:,i);
    vA_Des(:,i)=zeros(2,1);
    vA_Des_dot(:,i)=zeros(2,1);
    rAProj1(2*i-1:2*i,:)=rAProjO;
    SigmaProdD(i,1)=sigmaProdD;
end

%Control actions for the followers
end

function [uA_pot, minRAD]=potentialControl(XA,XD,rho_c_A,sigma_par,R_m_AD,R_bar_AD, R_u_AD, Rij,kADr,kADv, alphaADv)
global largeP

 rD=XD(1:2,:);
 vD=XD(3:4,:);
 rA=XA(1:2,1);
 vA=XA(3:4,1);
    
uADv=[0;0];
uADr=[0;0];

minRAD=Inf;
for j=1:length(rD(1,:))
    rADj=rD(1:2,j)-rA;
    R_ADj=norm(rADj);
    R_ADj_dot=(-rADj)'*vA/R_ADj;
    if R_ADj<=rho_c_A
        if R_ADj<R_u_AD
            if R_ADj<R_bar_AD
                sigma=1;     sigma_dot=0;
%                 ind=find(WA(i,:)==1);
%                 SigmaProdD(ind)=SigmaProdD(ind)*(1-sigma);
            elseif R_ADj>R_bar_AD && R_ADj<R_u_AD
                sigma=sigma_par(1)*R_ADj^3+sigma_par(2)*R_ADj^2+sigma_par(3)*R_ADj+sigma_par(4);
                %sigma_dot=(3*A_A_D*R_ADj^2+2*B_A_D*R_ADj+C_A_D)*R_ADj_dot;
            end
            %Rij=Rij0(j);
            if (R_ADj-R_m_AD)>1e-1
                nabla_ri_Vij=kADr*(rA-rD(1:2,j))/R_ADj/abs(R_ADj-R_m_AD)*((R_ADj-R_m_AD)^2-Rij^2)/((R_ADj-R_m_AD)^2+Rij^2);
            else
                nabla_ri_Vij=-kADr*(rA-rD(1:2,j))*largeP;
            end
            norm_dvA=norm(vA-vD(:,j));
            if norm_dvA>1e-16
            dv=kADv*(vA-vD(:,j))* norm_dvA^(alphaADv-1);
            else
                dv=[0,0]';
            end
            uADv=uADv-sigma*dv;
            uADr=uADr-sigma*nabla_ri_Vij;
        else
            sigma=0;   sigma_dot=0;
        end
%         Sigma(NO+j)=sigma;
%         Sigma_dot(NO+j)=sigma_dot;
%         F_ADj=(-rADj)/R_ADj;
%         F_AD(:,i) = F_AD(:,i) + sigma*F_ADj;
%         F_AD_dot(:,i) = F_AD_dot(:,i) + sigma_dot*F_ADj+sigma*(vA-R_ADj_dot*F_ADj)/R_ADj;
%         sigmaProdD=sigmaProdD*(1-sigma);
%         sigmaSumD=sigmaSumD+sigma;
    end
    if minRAD>R_ADj
        minRAD=R_ADj;
    end    
end
 uA_pot=uADv+uADr;
end
