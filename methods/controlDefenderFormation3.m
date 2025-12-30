function [ uD, arr_EDO, flagAttInSight] = controlDefenderFormation3(XD, indDef, assign, XD_des, XD_des_dot,uDFc_trans, WD, Rjj_bar, XA, NA, ND)
% ------------------------------
% Author: Vishnu S. Chipade
% Email: vishnuc@umich.edu
% ------------------------------

%This function finds control for defenders to maintain an open formation
%and move close to the attackers while the mouth of the formation faces the
%attackers

global obs rS rP nO aO bO kr0 kr1 kr2 kv1 kv2 RD_con v_maxD u_maxD u_maxA C_d rho_c_A options
global E_m_DO E_bar_DO E_u_DO E_v_DO R_m_DD R_bar_DD R_u_DD R_u_AD R_m_AD R_bar_AD
global A_D_O B_D_O C_D_O D_D_O A_D_D B_D_D C_D_D D_D_D
global A_bar_D_O B_bar_D_O C_bar_D_O D_bar_D_O
global A_A_D B_A_D C_A_D D_A_D
global sigmaProdD_dot_fun
global rSD_goal largeP
global kDOr kDDr kDOv kDOv2 alphaDOv kDDv alphaDDv kDFr kDFr2 kDFv alphaDFr alphaDFv kDRr kDRv
global Rjk00
global rho_Acon rho_sn
global kDFphi kDFphiv kDFphir kDFphid

na=length(assign);
nid=length(indDef);
indDef([assign;na+1:nid])=indDef;

XD=XD(:,indDef);
%XD_des=XD_des(:,indDef);
NO=obs.NO;
rCO2=obs.rCO2;

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


%Add the virutal defender centered at the CoM of the attackers
%XD=[XD,[rAcm;vAcm]];

rDFcm=XD(1:2,ND+1);
vDFcm=XD(3:4,ND+1);


%Find the control for CoM of the defenders by considering unicycle dynamcis
%for it and then tranform this control to the double integrator dynamcis
rDc_res=[0,0]'; %resultant of vector w.r.t formation center
vDc_res=[0,0]'; %resultant of vector w.r.t formation center
rho_DF0=zeros(ND,1);
for j=1:ND
    drDF=rDFcm-XD(1:2,j);
    rho_DF0(j)=norm(drDF);
    rDc_res=rDc_res+drDF;
    vDc_res=vDc_res+(vDFcm-XD(3:4,j));
end
rho_DF=max(rho_DF0);
drcm=rAcm-rDFcm;
dvcm=vAcm-vDFcm;
phiD_cm=atan2(rDc_res(2),rDc_res(1));
phiD_cm_dot=[rDc_res(1)*vDc_res(2)-rDc_res(2)*vDc_res(1)]/norm(rDc_res)^2;
phiD_cm_des=atan2(drcm(2),drcm(1));
phiD_cm_des_dot=[drcm(1)*dvcm(2)-drcm(2)*dvcm(1)]/norm(drcm)^2;
%phiD_cm_dot=-kDFphi*(phiD_cm-phiD_cm_des);
dphi=(phiD_cm-phiD_cm_des);
if dphi<-pi
    dphi=dphi+2*pi;
end
if dphi<10*pi/180
    flagAttInSight=1;
else
    flagAttInSight=0;
end
phiD_cm_ddot=-kDFphi*dphi-kDFphid*(phiD_cm_dot-phiD_cm_des_dot);

%vD_cm_dot=kDFphiv*(1/norm(rAcm-rDcm)*(rAcm-rDcm)'*(vAcm-vDcm)); %scalar acceleration (vD_cm=kDFphiv*norm(rAcm-rDcm))
% vD_cm_dot=kDFphir*(norm(rAcm-rDcm))-kDFphiv*norm(vAcm-vDcm);
% uD_cm=[cos(phiD_cm), -sin(phiD_cm); sin(phiD_cm),  cos(phiD_cm)]*[vD_cm_dot, norm(vDcm)*phiD_cm_dot]';
if norm(rDFcm-rAcm)<3*rho_sn
    if vAcm'*(rP-rDFcm)<0
        uDFcm0=-kDFphir*(rDFcm-rAcm)-kDFphiv*(vDFcm-vAcm);
    else
        uDFcm0=-kDFphir*(rDFcm-rAcm);
    end
else
    uDFcm0=-kDFphir*(rDFcm-rAcm);
end

%For obstacle avoidance
uDFcmOv=zeros(2,1);
uDFcmOr=uDFcmOv;
if(1)
    for k=1:NO
        n=nO(k);
        a=aO(k);
        b=bO(k);
        R_OD=norm(rDFcm-rCO2(1:2,k));
        EmDO=E_m_DO(1,k);
        EbarDO=E_bar_DO(1,k);
        EuDO=E_u_DO(1,k);
        EvDO=E_v_DO(1,k);
        E_DOk=(abs((rDFcm(1)-rCO2(1,k))/a))^(2*n)+(abs((rDFcm(2)-rCO2(2,k))/b))^(2*n)-1;
        %E_DOk_dot=2*n*[(abs((rDFcm(1)-rCO2(1,k))/a))^(2*n)/(rDFcm(1)-rCO2(1,k)), (abs((rDFcm(2)-rCO2(2,k))/b))^(2*n)/(rDFcm(2)-rCO2(2,k))]*vDcm;
        %arr_EDO(1,k)=E_DOk;
        if E_DOk<EuDO
            if E_DOk>EmDO &&  E_DOk<EbarDO
                sigma=1;
                %sigma_dot=0;
                sigma_bar=1;
            elseif  E_DOk>EbarDO &&  E_DOk<EuDO
                sigma=A_D_O(1,k)*E_DOk^3+B_D_O(1,k)*E_DOk^2+C_D_O(1,k)*E_DOk+D_D_O(1,k);
                %sigma_dot=(3*A_D_O(j,k)*E_DOk^2+2*B_D_O(j,k)*E_DOk+C_D_O(j,k))*E_DOk_dot;
                sigma_bar=1;
            end
            %Point along the radial direction as the initial condition
            beta=atan2(rDFcm(2)-rCO2(2,k),rDFcm(1)-rCO2(1,k));
            RProj=((EmDO+1)/(abs(cos(beta))^(2*n)/a^(2*n)+abs(sin(beta))^(2*n)/b^(2*n)))^(1/(2*n));
            rP0=[RProj*cos(beta);RProj*sin(beta)]+rCO2(1:2,k);
            f=@(r) [abs(r(1)-rCO2(1,k))^(2*n)/a^(2*n)+abs(r(2)-rCO2(2,k))^(2*n)/b^(2*n)-EmDO-1;...
                a^(2*n)*sign(r(2)-rCO2(2,k))*abs(r(2)-rCO2(2,k))^(2*n-1)*(rDFcm(1)-r(1))-b^(2*n)*sign(r(1)-rCO2(1,k))*abs(r(1)-rCO2(1,k))^(2*n-1)*(rDFcm(2)-r(2))];
            rDFcmProj(:,k)=fsolve(f,rP0,options);
            %Tangent at the projection point
            beta_barP=atan2(b^(2*n)*sign(rDFcmProj(1,k)-rCO2(1,k))*abs(rDFcmProj(1,k)-rCO2(1,k))^(2*n-1),-a^(2*n)*sign(rDFcmProj(2,k)-rCO2(2,k))*abs(rDFcmProj(2,k)-rCO2(2,k))^(2*n-1));
            rTP=[cos(beta_barP);sin(beta_barP)];
            if rTP'*vDFcm<0  %If the velocity of A is opposite to the tangent then reveres the tangent direction for projection
                rTP=-rTP;
            end
            vDFcmProj(:,k)=rTP'*vDFcm*rTP;
        elseif E_DOk>EuDO &&  E_DOk<EvDO
            sigma_bar=A_bar_D_O(1,k)*E_DOk^3+B_bar_D_O(1,k)*E_DOk^2+C_bar_D_O(1,k)*E_DOk+D_bar_D_O(1,k);
            sigma=0;
            %sigma_dot=0;
        else
            sigma=0;
            %sigma_dot=0;
            sigma_bar=0;
        end
        Sigma(k)=sigma;
        %Sigma_dot(k)=sigma_dot;
        if Sigma(k)~=0
            Rjk0=Rjk00(1);
            Rj_jk=norm(rDFcm-rDFcmProj(:,k));
            nabla_rj_VjP= kDOr*Sigma(k)*(rDFcm-rDFcmProj(:,k))/Rj_jk/abs(Rj_jk-(R_m_DD+rho_DF))*((Rj_jk-(R_m_DD+rho_DF))^2-Rjk0^2)/((Rj_jk-(R_m_DD+rho_DF))^2+Rjk0^2);
            dv=kDOv2*(vDFcm-vDFcmProj(:,k))*norm(vDFcm-vDFcmProj(:,k))^(alphaDOv-1);
            uDFcmOv=uDFcmOv-Sigma(k)*dv;
            uDFcmOr=uDFcmOr-Sigma(k)*nabla_rj_VjP;
        end
    end
end

uDFcm_trans=uDFcm0+uDFcmOr+uDFcmOv;
chi=0.7;
%Apply saturation to defenders' control action
infNorm_uDFcm_trans=max(abs(uDFcm_trans(1)),abs(uDFcm_trans(2)));
if infNorm_uDFcm_trans>chi*u_maxD(1)
    uDFcm_trans=uDFcm_trans*chi*u_maxD(1)/infNorm_uDFcm_trans;
end



%Since the entire formation is to be kept rigid, all defenders apply same
%control law
for j=1:ND
    uD_rot=phiD_cm_ddot*[-(XD(2,j)-rDFcm(2)),(XD(1,j)-rDFcm(1))]';
    %Apply saturation to defenders' control action
% infNorm_uD_rot=max(abs(uD_rot(1)),abs(uD_rot(2)));
% if infNorm_uD_rot>(1-chi)*u_maxD(1)
%     uD_rot=uD_rot*(1-chi)*u_maxD(1)/infNorm_uD_rot;
% end
    uD(:,j)=uDFcm_trans+uD_rot;
    %Find the distance to obstacles
    rD=XD(1:2,j);
    for k=1:NO
        n=nO(k);
        a=aO(k);
        b=bO(k);
        E_DOk=(abs((rD(1)-rCO2(1,k))/a))^(2*n)+(abs((rD(2)-rCO2(2,k))/b))^(2*n)-1;
        arr_EDO(j,k)=E_DOk;
    end
end
%the control for the virtual agent at the formation center
uD(:,ND+1)=uDFcm_trans;
%Apply saturation to defenders' control action
for j=1:ND
    %norm_uD=norm(uD(:,j));
    infNorm_uD=max(abs(uD(1,j)),abs(uD(2,j)));
    if infNorm_uD>u_maxD(j)
        uD(:,j)=uD(:,j)*u_maxD(j)/infNorm_uD;
    end
end

%Change the control according to the indices
uD(:,indDef)=uD;
end

