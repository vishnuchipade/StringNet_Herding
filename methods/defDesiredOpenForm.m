function [XD_des,XD_des_dot,phi_dot,phi_ddot,uDFc_trans,flagAttInSight]=defDesiredOpenForm(XDFc,RDF0,XA,phi,phi_dot,NA,ND)
% ------------------------------
% Author: Vishnu S. Chipade
% Email: vishnuc@umich.edu
% ------------------------------


global obs rP nO aO bO rho_sn u_maxD kDFphi kDFphid kDFphir kDFphiv
global E_m_DO E_bar_DO E_u_DO E_v_DO R_m_DD R_bar_DD R_u_DD R_u_AD R_m_AD R_bar_AD
global A_D_O B_D_O C_D_O D_D_O A_D_D B_D_D C_D_D D_D_D
global A_bar_D_O B_bar_D_O C_bar_D_O D_bar_D_O
global A_A_D B_A_D C_A_D D_A_D
global sigmaProdD_dot_fun
global rSD_goal largeP
global kDOr kDDr kDOv kDOv2 alphaDOv kDDv alphaDDv kDFr kDFr2 kDFv alphaDFr alphaDFv kDRr kDRv
global Rjk00
global options
global rho_Acon rhoD_safe
global R_bar_DO R_u_DO R_m_DDO
global C_d umdf_s1 umdf_s2

NO=obs.NO;
rCO2=obs.rCO2;
rVO=obs.rVO1;

chi=0.55;
phi_ddot_max=(0.9-chi)*u_maxD(1)/RDF0;


%Initial reference semicircular formation
rAcm=sum(XA(1:2,:),2)/NA;
vAcm=sum(XA(3:4,:),2)/NA;
rDFc=XDFc(1:2,1);
vDFc=XDFc(3:4,1);
drc=rDFc-rAcm;
dvc=vDFc-vAcm;

phi_des=atan2(-drc(2),-drc(1));
if phi_des<0
    phi_des=phi_des+2*pi;
end
phi_des_dot=[-drc(1)*dvc(2)+drc(2)*dvc(1)]/norm(drc)^2;
dphi=(phi-phi_des);
if dphi<-pi
    dphi=dphi+2*pi;
end
if dphi<abs(asin((RDF0)/norm(drc)))  %-rho_Acon
    flagAttInSight=1;
else
    flagAttInSight=0;
end
phi_ddot0=-kDFphi*dphi-kDFphid*(phi_dot-phi_des_dot);
phi_ddot=sign(phi_ddot0)*min(abs(phi_ddot0),phi_ddot_max);

%Translational part
% if norm(drc)<3*rho_sn
%     uDFc0=-kDFphir*(drc)-kDFphiv*(dvc);
% %     if vAcm'*(rP-rDFc)<0
% %         uDFc0=-kDFphir*(rDFc-rAcm)-kDFphiv*(vDFc-vAcm);
% %     else
% %         uDFc0=-kDFphir*(rDFc-rAcm);
% %     end
% else
%     uDFc0=-kDFphir*(rDFc-rAcm);
% end

uDFcr0=-kDFphir*(drc);
%uDFcv0=-kDFphiv*(dvc);

%For obstacle avoidance
uDFcOv=zeros(2,1);
uDFcOr=uDFcOv;
if(1)
    for k=1:NO
       R_DFcOk=Inf;
            %rVO{k}=[rVO{k},rVO{k}(:,1)];
            for kk=1:length(rVO{k}(1,:))-1
                dx=rVO{k}(1,kk+1)-rVO{k}(1,kk);
                dy=rVO{k}(2,kk+1)-rVO{k}(2,kk);
                if abs(dx)>1e-10
                    mVV=dy/dx;  %slope of the line
                    cVV=rVO{k}(2,kk)-mVV*rVO{k}(1,kk);
                    rDFcProjO_temp(1,1)=(mVV*rDFc(2)+rDFc(1)-mVV*cVV)/(1+mVV^2);
                    rDFcProjO_temp(2,1)=mVV*rDFcProjO_temp(1,1)+cVV;
                    lambdaDP=(rDFcProjO_temp(1,1)-rVO{k}(1,kk))/dx;
                else
                    rDFcProjO_temp(1,1)=rVO{k}(1,kk);
                    rDFcProjO_temp(2,1)=rDFc(2);
                    lambdaDP=(rDFcProjO_temp(2,1)-rVO{k}(2,kk))/dy;
                end
                if lambdaDP<0
                    rDFcProjO_temp(:,1)=rVO{k}(1:2,kk);
                elseif lambdaDP>1
                    rDFcProjO_temp(:,1)=rVO{k}(1:2,kk+1);
                end
                R_DFcOk_temp=norm(rDFc-rDFcProjO_temp);
                if R_DFcOk_temp<R_DFcOk
                    R_DFcOk=R_DFcOk_temp;
                    rDFcProjO(:,k)=rDFcProjO_temp;
                end
            end
            
            R_DFcOk=R_DFcOk-rhoD_safe;
            theta=atan2(rDFc(2)-rDFcProjO(2,k),rDFc(1)-rDFcProjO(1,k));
            rTP=[cos(theta+pi/2), sin(theta+pi/2)]';
            if rTP'*vDFc<0  %If the velocity of D is opposite to the tangent then reveres the tangent direction for projection
                rTP=-rTP;
            end
            vDFcProjO(:,k)=rTP'*vDFc*rTP;
           if R_DFcOk<R_u_DO(1,k)+RDF0
                if R_DFcOk<R_bar_DO(1,k)+RDF0  %&& E_AOk>EmAO
                    sigma=1;
                    sigma_dot=0;
                    sigma_bar=1;
                    FlagAttInObs(k)=1;
                elseif  R_DFcOk>R_bar_DO(1,k)+RDF0 &&  R_DFcOk<R_u_DO(1,k)+RDF0
                    sigma=A_D_O(1,k)*(R_DFcOk-RDF0)^3+B_D_O(1,k)*(R_DFcOk-RDF0)^2+C_D_O(1,k)*(R_DFcOk-RDF0)+D_D_O(1,k);
                   
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
            if (R_DFcOk-(R_m_DDO+RDF0))>1e-1
            nabla_rj_VjP= kDOr*Sigma(k)*(rDFc-rDFcProjO(:,k))/R_DFcOk/abs(R_DFcOk-(R_m_DDO+RDF0))*((R_DFcOk-(R_m_DDO+RDF0))^2-Rjk0^2)/((R_DFcOk-(R_m_DDO+RDF0))^2+Rjk0^2);
            else
            nabla_rj_VjP= - kDOr*Sigma(k)*(rDFc-rDFcProjO(:,k))*largeP;
            end
            dv=kDOv2*(vDFc-vDFcProjO(:,k))*norm(vDFc-vDFcProjO(:,k))^(alphaDOv-1);
            uDFcOv=uDFcOv-Sigma(k)*dv;
            uDFcOr=uDFcOr-Sigma(k)*nabla_rj_VjP;
        end
    end
end
uDFc_trans1=uDFcr0+ uDFcOr+uDFcOv;
norm_uDFc_trans1=norm(uDFc_trans1);
uDFc_trans2=-kDFphiv*(dvc)+C_d*vDFc*norm(vDFc);
norm_uDFc_trans2=norm(uDFc_trans2);

%Apply saturation to defenders' control action
uDFc_trans=min(umdf_s1,norm_uDFc_trans1)*uDFc_trans1/norm_uDFc_trans1+min(umdf_s2,norm_uDFc_trans2)*uDFc_trans2/norm_uDFc_trans2;

% infNorm_uDFc_trans=max(abs(uDFc_trans(1)),abs(uDFc_trans(2)));
% if infNorm_uDFc_trans>chi*u_maxD(1)
%     uDFc_trans=uDFc_trans*chi*u_maxD(1)/infNorm_uDFc_trans;
% end

%uDFc_trans=[0,0]';
%phi_ddot=0;
for j=1:ND
%     if flagNDeven~=1
%         phi_j=phi+2*pi*(j-1)/ND;
%     else
%         phi_j=phi+2*pi*(j-1)/ND+pi/ND;
%     end
     phi_j=phi+pi/2;%+(j-1)*pi/(ND-1);
    uD_rot=RDF0*phi_ddot*[-sin(phi_j),cos(phi_j)]'-RDF0*phi_dot^2*[cos(phi_j),sin(phi_j)]';    
    uD(:,j)=uDFc_trans+uD_rot;
    XD_des(1:2,j)=rDFc+(0.5*(ND-1)*RDF0-(j-1)*RDF0)*[cos(phi_j), sin(phi_j)]';
    XD_des(3:4,j)=vDFc+(0.5*(ND-1)*RDF0-(j-1)*RDF0)*phi_dot*[-sin(phi_j), cos(phi_j)]';
    XD_des_dot(3:4,j)=uD(:,j);
    XD_des_dot(1:2,j)=XD_des(3:4,j);
end
end