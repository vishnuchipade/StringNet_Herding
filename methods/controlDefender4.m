function [ uD, arr_EDO] = controlDefender4(XD,indDef, XD_des, XD_des_dot, XA, WD, Rjj_bar, ND )
% ------------------------------
% Author: Vishnu S. Chipade
% Email: vishnuc@umich.edu
% ------------------------------

%Flocking type control for defenders to herd the trapped attackers to the
%safe area
global obs rS nO aO bO kr0 kr1 kr2 kv1 kv2 RD_con v_maxD u_maxD u_maxA C_d rho_c_A options
global E_m_DcO E_bar_DcO E_u_DcO E_v_DcO R_m_DD R_bar_DD R_u_DD R_u_AD R_m_AD R_bar_AD
global A_Dc_O B_Dc_O C_Dc_O D_Dc_O A_D_D B_D_D C_D_D D_D_D
global A_bar_Dc_O B_bar_Dc_O C_bar_Dc_O D_bar_Dc_O
global A_A_D B_A_D C_A_D D_A_D
global sigmaProdD_dot_fun
global rSD_goal largeP
global kDOr kDDr kDOv1 alphaDOv kDDv alphaDDv
global Rjk00 rho_sn Rjj0

XD=XD(:,indDef);
Rji0=0.6*R_bar_AD*ones(1,ND);

%indAtt=pairDA(indDef);
%rACon=XA(1:2,indAtt);
%vACon=XA(3:4,indAtt);

rDcm=sum(XD(1:2,:),2)/ND;
vDcm=sum(XD(3:4,:),2)/ND;

rA=XA(1:2,:);
vA=XA(3:4,:);
NA=size(XA,2);
NO=obs.NO;
rCO2=obs.rCO2;


%Center of mass of defenders
uDcmOv=zeros(2,1);
uDcmOr=uDcmOv;
%Point along the radial direction as the initial condition
%remember all the j's are changed to 1 for the CoM
for k=1:NO
    n=nO(k);
    a=aO(k);
    b=bO(k);
    R_OD=norm(rDcm-rCO2(1:2,k));
    EmDO=E_m_DcO(1,k);
    EbarDO=E_bar_DcO(1,k);
    EuDO=E_u_DcO(1,k);
    EvDO=E_v_DcO(1,k);
    E_DOk=(abs((rDcm(1)-rCO2(1,k))/a))^(2*n)+(abs((rDcm(2)-rCO2(2,k))/b))^(2*n)-1;
    %E_DOk_dot=2*n*[(abs((rDcm(1)-rCO2(1,k))/a))^(2*n)/(rDcm(1)-rCO2(1,k)), (abs((rDcm(2)-rCO2(2,k))/b))^(2*n)/(rDcm(2)-rCO2(2,k))]*vDcm;
    %arr_EDO(1,k)=E_DOk;
    if E_DOk<EuDO
        if E_DOk>EmDO &&  E_DOk<EbarDO
            sigma=1;
            %sigma_dot=0;
            sigma_bar=1;
        elseif  E_DOk>EbarDO &&  E_DOk<EuDO
            sigma=A_Dc_O(1,k)*E_DOk^3+B_Dc_O(1,k)*E_DOk^2+C_Dc_O(1,k)*E_DOk+D_Dc_O(1,k);
            %sigma_dot=(3*A_D_O(j,k)*E_DOk^2+2*B_D_O(j,k)*E_DOk+C_D_O(j,k))*E_DOk_dot;
            sigma_bar=1;
        else
            sigma=0;
        end
        %Point along the radial direction as the initial condition
        beta=atan2(rDcm(2)-rCO2(2,k),rDcm(1)-rCO2(1,k));
        RProj=((EmDO+1)/(abs(cos(beta))^(2*n)/a^(2*n)+abs(sin(beta))^(2*n)/b^(2*n)))^(1/(2*n));
        rP0=[RProj*cos(beta);RProj*sin(beta)]+rCO2(1:2,k);
        f=@(r) [abs(r(1)-rCO2(1,k))^(2*n)/a^(2*n)+abs(r(2)-rCO2(2,k))^(2*n)/b^(2*n)-EmDO-1;...
            a^(2*n)*sign(r(2)-rCO2(2,k))*abs(r(2)-rCO2(2,k))^(2*n-1)*(rDcm(1)-r(1))-b^(2*n)*sign(r(1)-rCO2(1,k))*abs(r(1)-rCO2(1,k))^(2*n-1)*(rDcm(2)-r(2))];
        rDcmProj(:,k)=fsolve(f,rP0,options);
        %Tangent at the projection point
        beta_barP=atan2(b^(2*n)*sign(rDcmProj(1,k)-rCO2(1,k))*abs(rDcmProj(1,k)-rCO2(1,k))^(2*n-1),-a^(2*n)*sign(rDcmProj(2,k)-rCO2(2,k))*abs(rDcmProj(2,k)-rCO2(2,k))^(2*n-1));
        rTP=[cos(beta_barP);sin(beta_barP)];
        if rTP'*vDcm<0  %If the velocity of A is opposite to the tangent then reveres the tangent direction for projection
            rTP=-rTP;
        end
        vDcmProj(:,k)=rTP'*vDcm*rTP;
    elseif E_DOk>EuDO &&  E_DOk<EvDO
        sigma_bar=A_bar_Dc_O(1,k)*E_DOk^3+B_bar_Dc_O(1,k)*E_DOk^2+C_bar_Dc_O(1,k)*E_DOk+D_bar_Dc_O(1,k);
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
        Rjk0=1.5*rho_sn;
        Rj_jk=norm(rDcm-rDcmProj(:,k));
        nabla_rj_VjP= kDOr*Sigma(k)*(rDcm-rDcmProj(:,k))/Rj_jk/abs(Rj_jk-rho_sn)*((Rj_jk-rho_sn)^2-Rjk0^2)/((Rj_jk-rho_sn)^2+Rjk0^2);
        dv=kDOv1*(vDcm-vDcmProj(:,k))*norm(vDcm-vDcmProj(:,k))^(alphaDOv-1);
        uDcmOv=uDcmOv-Sigma(k)*dv;
        uDcmOr=uDcmOr-Sigma(k)*nabla_rj_VjP;
    end
end
uD(:,ND+1)=[0,0]';

for j=1:ND
    rD=XD(1:2,j);
    vD=XD(3:4,j);
    Rjj=Rjj_bar(j,ND+1);
    rDjDc=rDcm-rD;
    R_DjDc=norm(rDcm-rD);
    if abs(R_DjDc-R_m_DD)>1e-10
        nabla_rj_Vjc=kDDr*(-rDjDc)/R_DjDc/abs(R_DjDc-R_m_DD)*((R_DjDc-R_m_DD)^2-Rjj^2)/((R_DjDc-R_m_DD)^2+Rjj^2);
    else
        nabla_rj_Vjc=kDDr*(-rDjDc)*largeP;
    end
    uD(:,j)=-0.1*vDcm-0.0005*(rDcm-rS)*norm(rDcm-rS)^(0.5-1)+uDcmOv+uDcmOr;%-nabla_rj_Vjc;  %0.012
    %Find the distance to obstacles
    %rD=XD(1:2,j);
    for k=1:NO
        n=nO(k);
        a=aO(k);
        b=bO(k);
        E_DOk=(abs((rD(1)-rCO2(1,k))/a))^(2*n)+(abs((rD(2)-rCO2(2,k))/b))^(2*n)-1;
        arr_EDO(j,k)=E_DOk;        
    end
end


%Apply saturation to defenders' control action
%This should be upper bounded by the attackers' acceleration
for j=1:ND
    %norm_uD=norm(uD(:,j));
    infNorm_uD=max(abs(uD(1,j)),abs(uD(2,j)));
    if infNorm_uD>0.95*u_maxA(1)
        uD(:,j)=uD(:,j)*0.95*u_maxA(1)/infNorm_uD;
    end
end

uD(:,indDef)=uD;
end
