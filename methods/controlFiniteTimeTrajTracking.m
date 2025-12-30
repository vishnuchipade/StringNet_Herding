function [ uD, arr_RDO] = controlFiniteTimeTrajTracking(XD,indDef, XD_des, XD_des_dot,uDFc_trans, XA, ND, flagAvoidACon )
% ------------------------------
% Author: Vishnu S. Chipade
% Email: vishnuc@umich.edu
% ------------------------------

%Defenders track their desired positions while avoiding obstacles and the
%connectivity region of the attackers
global obs C_d 
global R_m_DD R_m_DO  R_bar_DD R_u_DD R_bar_AD
global R_u_DO R_bar_DO
global A_D_O B_D_O C_D_O D_D_O A_D_D B_D_D C_D_D D_D_D
global rSD_goal largeP
global kDOr kDDr kDOv1 kDOv2 alphaDOv kDDv alphaDDv kDFr2 kDFv alphaDFr alphaDFv
global Rjk00 Rjj0
global rho_Acon
global rhoD_safe
global umd1 umd2

% na=length(assign);
% nid=length(indDef);
% indDef([assign;na+1:nid])=indDef;

XD=XD(:,indDef);
NA=size(XA,2);
rCO2=obs.rCO2;
NO=size(rCO2,2);
rVO=obs.rVO;

Rji0=0.6*R_bar_AD*ones(1,ND);

%indAtt=pairDA(indDef);
%rACon=XA(1:2,indAtt);
%vACon=XA(3:4,indAtt);

rDcm=sum(XD(1:2,:),2)/ND;
vDcm=sum(XD(3:4,:),2)/ND;

rA=XA(1:2,:);
vA=XA(3:4,:);

rAcm=sum(XA(1:2,:),2)/NA;
vAcm=sum(XA(3:4,:),2)/NA;

for j=1:ND
    rD=XD(1:2,j);
    rD_des=XD_des(1:2,j);
    vD_des=vAcm;
    vD=XD(3:4,j);
    
    %For converging to the formation
    dr0= rD-rD_des + 1/(kDFr2*(2-alphaDFv))*(vD-vD_des)*norm(vD-vD_des)^(1-alphaDFv);
    if norm(dr0)>1e-5
        dr=-kDFr2*(dr0)*norm(dr0)^(alphaDFr-1);
    else
        dr=[0,0]';
    end
    if norm(vD-vD_des)>1e-5
        dv=-kDFr2*(vD-vD_des)*norm(vD-vD_des)^(alphaDFv-1);
    else
        dv=[0,0]';
    end
    
    %check for nearby defenders
    uDDv=zeros(2,1);
    uDDr=zeros(2,1);
    if(1)
        count=0;
        for l=1:ND
            if l~=j
                count=count+1;
                rDjDl=XD(1:2,l)-rD;
                R_DjDl=norm(rDjDl);
                R_DjDl_dot=(-rDjDl)'*vD/R_DjDl;
                if R_DjDl<R_u_DD
                    if  R_DjDl<R_bar_DD %R_DjDl>R_m_DD &&
                        sigma=1;
                        sigma_dot=0;
                    elseif R_DjDl>R_bar_DD && R_DjDl<R_u_DD
                        sigma=A_D_D*R_DjDl^3+B_D_D*R_DjDl^2+C_D_D*R_DjDl+D_D_D;
                        sigma_dot=(3*A_D_D*R_DjDl^2+B_D_D*R_DjDl+C_D_D)*R_DjDl_dot;
                    end
                    Rjj=Rjj0(j);
                    nabla_rj_Vjj=kDDr*(-rDjDl)/R_DjDl/abs(R_DjDl-R_m_DD)*((R_DjDl-R_m_DD)^2-Rjj^2)/((R_DjDl-R_m_DD)^2+Rjj^2);
                    dv=kDDv*(vD-XD(3:4,l))*norm(vD-XD(3:4,l))^(alphaDDv-1);
                    uDDv=uDDv-sigma*dv;
                    uDDr=uDDr-sigma*nabla_rj_Vjj;
                else
                    sigma=0;
                    sigma_dot=0;
                end
                Sigma(count)=sigma;
                Sigma_dot(count)=sigma_dot;
                
                %                 sigmaProd=sigmaProd*(1-sigma);
                %                 sigmaSum=sigmaSum+sigma;
            end
        end
        %sigmaSum
    end
    
    %Avoid the boundary of the connectivity region of the attackers
    nabla_rj_VjP=[0,0]';
    dvP=[0,0]';
    if ( flagAvoidACon )
        thetaDAcm=atan2(rD(2)-rAcm(2),rD(1)-rAcm(1));
        rDProj=rAcm+rho_Acon*[cos(thetaDAcm),sin(thetaDAcm)]';
        rTP=[cos(thetaDAcm+pi/2),sin(thetaDAcm+pi/2)]';
        if rTP'*vD<0  %If the velocity of A is opposite to the tangent then reveres the tangent direction for projection
            rTP=-rTP;
        end
        vDProj=rTP'*vD*rTP+vAcm;
        Rjk0=Rjk00(1);
        Rj_jk=norm(rD-rDProj);
        if (Rj_jk-R_m_DD)>1e-1
            nabla_rj_VjP= kDOr*(rD-rDProj)/Rj_jk/abs(Rj_jk-R_m_DD)*((Rj_jk-R_m_DD)^2-Rjk0^2)/((Rj_jk-R_m_DD)^2+Rjk0^2);
        else
            nabla_rj_VjP= -kDOr*(rD-rDProj)*largeP;
        end
        if norm(vD-vDProj)~=0
            dvP=-kDOv2*(vD-vDProj)*norm(vD-vDProj)^(alphaDOv-1);
        else
            dvP=[0,0]';
        end
    end
    
    uDOv=zeros(2,1);
uDOr=uDOv;
if(1)
    for k=1:NO
        R_DOk=Inf;
            %rVO{k}=[rVO{k},rVO{k}(:,1)];
            for kk=1:length(rVO{k}(1,:))-1
                dx=rVO{k}(1,kk+1)-rVO{k}(1,kk);
                dy=rVO{k}(2,kk+1)-rVO{k}(2,kk);
                if abs(dx)>1e-10
                    mVV=dy/dx;  %slope of the line
                    cVV=rVO{k}(2,kk)-mVV*rVO{k}(1,kk);
                    rDProjO_temp(1,1)=(mVV*rD(2)+rD(1)-mVV*cVV)/(1+mVV^2);
                    rDProjO_temp(2,1)=mVV*rDProjO_temp(1,1)+cVV;
                    lambdaDP=(rDProjO_temp(1,1)-rVO{k}(1,kk))/dx;
                else
                    rDProjO_temp(1,1)=rVO{k}(1,kk);
                    rDProjO_temp(2,1)=rD(2);
                    lambdaDP=(rDProjO_temp(2,1)-rVO{k}(2,kk))/dy;
                end
                if lambdaDP<0
                    rDProjO_temp(:,1)=rVO{k}(1:2,kk);
                elseif lambdaDP>1
                    rDProjO_temp(:,1)=rVO{k}(1:2,kk+1);
                end
                R_DOk_temp=norm(rD-rDProjO_temp);
                if R_DOk_temp<R_DOk
                    R_DOk=R_DOk_temp;
                    rDProjO(:,k)=rDProjO_temp;
                end
            end
            
            R_DOk=R_DOk-rhoD_safe;
            arr_RDO(j,k)=R_DOk;
            theta=atan2(rD(2)-rDProjO(2,k),rD(1)-rDProjO(1,k));
            rTP=[cos(theta+pi/2), sin(theta+pi/2)]';
            if rTP'*vD<0  %If the velocity of D is opposite to the tangent then reveres the tangent direction for projection
                rTP=-rTP;
            end
            vDProjO(:,k)=rTP'*vD*rTP;
            if R_DOk<R_u_DO(1,k)
                if R_DOk<R_bar_DO(1,k)  %&& E_AOk>EmAO
                    sigma=1;
                    sigma_dot=0;
                    sigma_bar=1;
                    FlagAttInObs(k)=1;
                elseif  R_DOk>R_bar_DO(1,k) &&  R_DOk<R_u_DO(1,k)
                    sigma=A_D_O(1,k)*(R_DOk)^3+B_D_O(1,k)*(R_DOk)^2+C_D_O(1,k)*(R_DOk)+D_D_O(1,k);
                   
                    sigma_bar=1;
                end
            else
                sigma=0;
                sigma_dot=0;
            end
        Sigma(k)=sigma;
        %Sigma_dot(k)=sigma_dot;
        if Sigma(k)~=0
            Rjk0=Rjk00(1);
            if R_DOk-(R_m_DD)>1e-1
            nabla_rj_VjPO= kDOr*Sigma(k)*(rD-rDProjO(:,k))/R_DOk/abs(R_DOk-(R_m_DD))*((R_DOk-(R_m_DD))^2-Rjk0^2)/((R_DOk-(R_m_DD))^2+Rjk0^2);
            else
            nabla_rj_VjPO= - kDOr*Sigma(k)*(rD-rDProjO(:,k))*largeP;
            end
            dvPO=kDOv2*(vD-vDProjO(:,k))*norm(vD-vDProjO(:,k))^(alphaDOv-1);
            uDOv=uDOv-Sigma(k)*dvPO;
            uDOr=uDOr-Sigma(k)*nabla_rj_VjPO;
        end
    end
    end
    
    
    
    %Apply saturation to defenders' control action
    uD1=dr+uDDv+uDDr+dvP-nabla_rj_VjP+uDOv+uDOr;
    norm_uD1=norm(uD1);
    uD2=dv+C_d*XD(3:4,j)*norm(XD(3:4,j));
    norm_uD2=norm(uD2);
    %Control Action for the defender
    if norm_uD1>1e-15
        uD1_hat=uD1/norm_uD1;
    else
        uD1_hat=[0,0]';
    end
    if norm_uD2>1e-15
        uD2_hat=uD2/norm_uD2;
    else
        uD2_hat=[0,0]';
    end
    uD(:,j)=XD_des_dot(3:4,j)+min(umd1,norm_uD1)*uD1_hat+min(umd2,norm_uD2)*uD2_hat;
end

%the control for the virtual agent at the formation center
uD(:,ND+1)=uDFc_trans;

%Change the control according to the indices
uD(:,indDef)=uD;
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
