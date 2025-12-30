function [ uA, uA0, R_AO_min, R_AAProjS_min, vA_Des, vA_Des_dot, rAProj1, FlagAttInObs, BetaAv0, SpeedA0,mAv0,cAv0, SigmaProdD,F_A, F_A_dot] = controlAttacker4( XA, XA_goal, XA_goal_dot,A_lead, FlagAttInObs,flagEnclose, flagHerd,BetaAv0,SpeedA0,mAv0,cAv0, XD, WA, WDString, NA, ND )
% ------------------------------
% Author: Vishnu S. Chipade
% Email: vishnuc@umich.edu
% ------------------------------

%Control action for a swarm of attackers, in some configuration given by
%their communication graph and inter-agent distance potential,
%to reach protected area while avoiding obstacles and defenders

%Obstacle avoidance using beta-agent idea rather than vector fields for the
%flock
global obs rP C_d rhoAD_safe rho_safe
global u_maxA v_maxA rho_c_A
global R_m_AO R_bar_AO R_u_AO R_m_AD R_bar_AD R_u_AD R_m_AA R_bar_AA R_u_AA
global A_A_O B_A_O C_A_O D_A_O A_A_D B_A_D C_A_D D_A_D A_A_A B_A_A C_A_A D_A_A
global sigmaProdA_dot_fun
global Rii_bar Rij0 Rik00 options
global kd largeP rho_Acon
global kADr kAFr kAOr kAOr2 kAPr kAPv kAOv2 kAFv alphaAFv kAOv alphaAOv kADv alphaADv
global NA_sep

global R_m_AD2 R_bar_AD2 R_u_AD2
global A_A_D2 B_A_D2 C_A_D2 D_A_D2

tol=0.5;

NO=obs.NO;
rCO2=obs.rCO2;
rVO=obs.rVO1;  %to have the first vertex repeated
if ND>0
    rD=XD(1:2,1:ND);
    vD=XD(3:4,1:ND);
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
            
            R_AOk=R_AOk-rho_safe;
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
                
                if (R_AOk-R_m_AA)>tol
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
    if(1)
        uAFv=zeros(2,1);
        uAFr=uAFv;
        for ii=find(WA(i,:)==1)
            Rii0=Rii_bar(i,ii)-R_m_AA;  %k
            Ri_ii=norm(rA-XA(1:2,ii));
            if Ri_ii<R_u_AA
                if Ri_ii<R_bar_AA  %&& E_AOk>EmAO
                    sigma=1;
                    sigma_dot=0;
                    sigma_bar=1;
                    FlagAttInObs(i,k)=1;
                elseif  Ri_ii>R_bar_AA &&  Ri_ii<R_u_AA
                    sigma=A_A_A*Ri_ii^3+B_A_A*Ri_ii^2+C_A_A*Ri_ii+D_A_A;
                    %sigma_dot=(3*A_A_O(i,k)*R_AOk^2+2*B_A_O(i,k)*R_AOk+C_A_O(i,k))*R_AOk_dot;
                    sigma_bar=1;
                end
            else
                sigma=0;
                sigma_dot=0;
                sigma_bar=0;
                FlagAttInObs(i,k)=0;
            end
            if Ri_ii-R_m_AA>tol
                nabla_ri_Vii=kAFr*WA(i,ii)*(rA-XA(1:2,ii))/Ri_ii/abs(Ri_ii-R_m_AA)*((Ri_ii-R_m_AA)^2-Rii0^2)/((Ri_ii-R_m_AA)^2+Rii0^2);
            else
                nabla_ri_Vii=-kAFr*WA(i,ii)*(rA-XA(1:2,ii))*largeP;
            end
            if norm(vA-XA(3:4,ii))>1e-16
                dv=kAFv*(vA-XA(3:4,ii))*norm(vA-XA(3:4,ii))^(alphaAFv-1);
            else
                dv=[0,0]';
            end
            uAFv=uAFv-sigma*dv;
            uAFr=uAFr-sigma*nabla_ri_Vii;
        end
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
        if i>NA-NA_sep
            R_m=15*R_m_AD;
        R_underbar=R_m+20;
        R_bar=R_m+25;  
            [uAD_pot, minRAD]=potentialControl(XA(:,i),XD(:,1:ND),2*rho_c_A,sigma_parameters(R_underbar,R_bar),R_m,R_underbar,R_bar, R_bar+10,kADr,kADv, alphaADv);
        else
            [uAD_pot, minRAD]=potentialControl(XA(:,i),XD(:,1:ND),rho_c_A,[A_A_D,B_A_D, C_A_D, D_A_D],R_m_AD,R_bar_AD, R_u_AD, Rij0(1),kADr,kADv, alphaADv);
        end
        %potential control correspnding to the strings.
         R_underbar = R_bar_AD;
         R_bar = R_u_AD;
         R_m = R_m_AD;
         if flagHerd == 1  % for MAir experiments if already hearding only apply action corresponding the the center of mass of the attackers
             rA_temp = rAcm;
             vA_temp = vAcm;
             R_underbar = R_bar_AD + rho_Acon;
             R_bar = R_u_AD + rho_Acon;
             R_m = 1*R_m_AD + rho_Acon;
             distMin=Inf;
             for j=1:ND
                 dist = norm(rAcm-XD(1:2,j));
                 if dist<distMin
                     XDClosest = XD(1:2,i);
                     distMin=dist;
                 end
             end
             uAD_pot = potentialControl([rAcm;vAcm],XD(:,1:ND),rho_c_A,sigma_parameters(R_underbar,R_bar),R_m,R_underbar,R_bar, R_bar+10,kADr,kADv, alphaADv);
         else
             rA_temp=rA;
             vA_temp=vA;
         end
        for j=1:ND                    
            %Find the projections of the attacker on the string
            %attached to defender j and find control corresponding to it
            uAAProj=zeros(2,1);
            for jj=find(WDString_temp(j,:)==1)
                countAPS = countAPS +1;
                [rAProjS(:,countAPS), lambdaAP]=projectionOnLine(rA_temp,XD(1:2,jj),XD(1:2,j));
                rTDD=XD(1:2,j)-XD(1:2,jj);
                rTDD=rTDD/norm(rTDD);
                if rTDD'*vA_temp<0
                    rTDD=-rTDD;
                end
                %velocity projection due to Attacker
                vAProjS0(:,countAPS)=rTDD'*vA_temp*rTDD;
                %velocity because of the motion of the defenders
                vAProjS1(:,countAPS)=XD(3:4,j)+(XD(3:4,jj)-XD(3:4,j))*norm(rAProjS(1:2,countAPS)-XD(1:2,j))/norm(XD(1:2,jj)-XD(1:2,j));
                vAProjS(:,countAPS)=vAProjS0(:,countAPS)+vAProjS1(:,countAPS);                %Remove the strings from calculations
                WDString_temp(j,jj)=0;
                WDString_temp(jj,j)=0;
            end
        end
        
        if flagHerd ==1
            if countAPS >2
            for js = 1:countAPS 
                dist(js)= norm(rAcm-rAProjS(1:2,js));
            end
            [values,ind]=sort(dist);
            rAProjS = rAProjS(:,ind(1:2)); %choose only the closest two strings
            vAProjS = vAProjS(:,ind(1:2));
            countAPS = 2;
            end
        end
        
        for js = 1:countAPS           
                uAAProj = uAAProj + potentialControl([rA_temp;vA_temp],[rAProjS(1:2,js);vAProjS(:,js)],2*rho_c_A,sigma_parameters(R_underbar,R_bar),R_m,R_underbar,R_bar, R_bar+10,kADr,kADv, alphaADv);
        end
    end
    
     
    %Stay within the connectivity region of the attackers
    if (1)
        if norm(rA-rAcm)<1.2*rho_Acon
            thetaAAcm=atan2(rA(2)-rAcm(2),rA(1)-rAcm(1));
            rAProjC=rAcm+rho_Acon*[cos(thetaAAcm),sin(thetaAAcm)]';
            rTP=[cos(thetaAAcm+pi/2),sin(thetaAAcm+pi/2)]';
            if rTP'*vA<0  %If the velocity of A is opposite to the tangent then reveres the tangent direction for projection
                rTP=-rTP;
            end
            vAProjC=rTP'*vA*rTP+vAcm;
            Rik0=Rik00(1);
            Ri_ik=norm(rA-rAProjC);
            if Ri_ik-R_m_AA>tol
                nabla_ri_ViP= kAOr2*(rA-rAProjC)/Ri_ik/abs(Ri_ik-R_m_AA)*((Ri_ik-R_m_AA)^2-Rik0^2)/((Ri_ik-R_m_AA)^2+Rik0^2);
            else
                nabla_ri_ViP= - kAOr2*(rA-rAProjC)*largeP;
            end
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
        if i==1
            duA_goal=XA_goal_dot(3:4,i)-(kAPr*(rA-rA_goal)+kAPv*(vA-vA_goal));  %sigmaProdD*sigmaProdDS*
        else
            duA_goal=XA_goal_dot(3:4,i)-(10*kAPr*(rA-rA_goal)+10*kAPv*(vA-vA_goal));
        end
        norm_duA_goal=norm(duA_goal);
        if (1)% A_lead(i)==2
            if norm_duA_goal>1e-10 && flagHerd~=1
                uA(:,i)=duA_goal+uAFv+uAFr+uAOv+uAOr+uAD_pot+uAAProj+C_d*norm(vA)*vA;% +dvP -nabla_ri_ViP; %+uADv
            else
                uA(:,i)=uAFv+uAFr+uAOv+uAOr+uAD_pot+uAAProj+C_d*norm(vA)*vA;% +uAFv+uAFr+uAFv+uAFr++dvP -nabla_ri_ViP; %+uADv
            end
        elseif (0)
            if norm_duA_goal>1e-10
                uA(:,i)=min(u_maxA(i),norm_duA_goal)*(duA_goal/norm_duA_goal)+uAFv+uAFr+uAOv+uAOr+uAD_pot+uAAProj+C_d*norm(vA)*vA;% +dvP -nabla_ri_ViP; %+uADv
            else
                uA(:,i)=uAFv+uAFr+uAOv+uAOr+uAD_pot+uAAProj+C_d*norm(vA)*vA;% +uAFv+uAFr+uAFv+uAFr++dvP -nabla_ri_ViP; %+uADv
            end
        else
            if norm_duA_goal>1e-10
                uA(:,i)=min(.8*u_maxA(i),norm_duA_goal)*(duA_goal/norm_duA_goal)+uAFv+uAFr+uAOv+uAOr+uAD_pot+uAAProj+C_d*norm(vA)*vA;% +dvP -nabla_ri_ViP; %+uADv
            else
                uA(:,i)=uAFv+uAFr+uAOv+uAOr+uAD_pot+uAAProj+C_d*norm(vA)*vA;% +uAFv+uAFr+uAFv+uAFr++dvP -nabla_ri_ViP; %+uADv
            end
        end
        uA0(:,i)=uA(:,i);
        % uA(:,i)=-sigmaProdD*(0.005*(rA-rA_goal)+0.5*(vA-vA_goal));
    else
        uA(:,i)=uA00(:,i)+uADv+uADr+uAP0(:,i);
    end
    if  1 && flagHerd~=1 && i==1  %flagEnclose~=1
        R_m=1.5*R_m_AD;
        R_underbar=R_m+2;
        R_bar=R_m+5;        
    [uAD_pot2, minRAD]=potentialControl(XA(:,i),XD(:,[2:ND]),2*rho_c_A,sigma_parameters(R_underbar,R_bar),R_m,R_underbar,R_bar, R_bar+10,kADr,kADv, alphaADv);

    uA(:,i)=uA(:,i)+uAD_pot2;
    end
    %infNorm_uA=max(abs(uA(1,i)),abs(uA(2,i)));
    
    
    %apply saturation
    norm_uA=norm(uA(:,i));
    
    if flagEnclose==1
        uMaxA = 0.9*u_maxA(i);
    else
        uMaxA = u_maxA(i);
    end
    if norm_uA>uMaxA %infNorm_uA>u_maxA(i)
        %uA(:,i)=uA(:,i)*u_maxA(i)/infNorm_uA;
        uA(:,i)=uA(:,i)*uMaxA/norm_uA;
    end
    if isnan(uA(1,i))
        1
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
%% Supporting functions
%Sigma parameters for smooth treansition of cotnrol action
    function sigma_params=sigma_parameters(R_underbar,R_bar)
        dR_cube=(R_bar-R_underbar)^3;
        sigma_params(1,1)=2/dR_cube;
        sigma_params(1,2)=-3*(R_bar+R_underbar)/dR_cube;
        sigma_params(1,3)=6*R_bar*R_underbar/dR_cube;
        sigma_params(1,4)=R_bar^2*(R_bar-3*R_underbar)/dR_cube;
    end

%Control for Collision avoidance of agent 1 from agents 2 using
%potential function approach
    function [u_pot_12, minDist_12]=potentialControl(X1,X2,rho_sens_1,sigma_params_12,R_m_12,R_underbar_12, R_bar_12, R_tilde_12,kr_12,kv_12, alphav_12)
        
        r2=X2(1:2,:);
        v2=X2(3:4,:);
        r1=X1(1:2,1);
        v1=X1(3:4,1);
        
        uv_12=[0;0];
        ur_12=[0;0];
        
        minDist_12=Inf;
        tol=1e-5;
        tol=0.1;
        largeNum=-1/tol*((tol)^2-R_tilde_12^2)/((tol)^2+R_tilde_12^2);
        for j=1:length(r2(1,:))
            r12j=r2(1:2,j)-r1;
            R_12j=norm(r12j);
            %R_ADj_dot=(-rADj)'*vA/R_ADj;
            if R_12j<=rho_sens_1
                if R_12j<R_bar_12
                    if R_12j<R_underbar_12
                        sigma=1;
                    elseif R_12j>R_underbar_12 && R_12j<R_bar_12
                        sigma=sigma_params_12(1)*R_12j^3+sigma_params_12(2)*R_12j^2+sigma_params_12(3)*R_12j+sigma_params_12(4);
                    end
                    if (R_12j-R_m_12)>tol
                        nabla_r1_V12j=kr_12*(r1-r2(1:2,j))/R_12j/abs(R_12j-R_m_12)*((R_12j-R_m_12)^2-R_tilde_12^2)/((R_12j-R_m_12)^2+R_tilde_12^2);
                    else
                        nabla_r1_V12j=-kr_12*(r1-r2(1:2,j))/R_12j*largeNum;
                    end
                    norm_dv12=norm(v1-v2(:,j));
                    if norm_dv12>1e-16 && ((v1-v2(:,j))'*(r1-r2(:,j))<0)
                        dv=kv_12*(v1-v2(:,j))* norm_dv12^(alphav_12-1);
                    else
                        dv=[0,0]';
                    end
                    uv_12=uv_12-sigma*dv;
                    ur_12=ur_12-sigma*nabla_r1_V12j;
                end
            end
            if minDist_12>R_12j
                minDist_12=R_12j;
            end
        end
        u_pot_12=uv_12+ur_12;
    end

    function [rProj,lambdaP]=projectionOnLine(r,r1,r2)
        dx=r2(1)-r1(1);
        dy=r2(2)-r1(2);
        mVV=dy/dx;  %slope of line joining the defenders
        cVV=r1(2)-mVV*r1(1);
        rProj(1,1)=(mVV*r(2)+r(1)-mVV*cVV)/(1+mVV^2);
        rProj(2,1)=mVV*rProj(1)+cVV;
        if mVV<1e16
            lambdaP=(rProj(1)-r1(1))/dx;
        else
            lambdaP=(rProj(2)-r1(2))/dy;
        end
        if lambdaP<0
            rProj=r1;
        elseif lambdaP>1
            rProj=r2;
        end
    end