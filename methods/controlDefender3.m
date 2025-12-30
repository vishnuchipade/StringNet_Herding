function [ uD, arr_EDO,FlagDefInObs,BetaDv0,mDv0,cDv0 ] = controlDefender3(XD,indDef, XD_des, XD_des_dot,FlagDefInObs,BetaDv0,mDv0,cDv0,XA,pairDA, refTraj, ND )
% ------------------------------
% Author: Vishnu S. Chipade
% Email: vishnuc@umich.edu
% ------------------------------

%Flocking type control for defenders
global rVO rCO2 rS nO aO bO kr0 kr1 kr2 kv1 kv2 RD_con v_maxD u_maxD C_d rho_c_A options
global E_m_DO E_bar_DO E_u_DO E_v_DO R_m_DD R_bar_DD R_u_DD R_u_AD R_m_AD R_bar_AD
global A_D_O B_D_O C_D_O D_D_O A_D_D B_D_D C_D_D D_D_D
global A_bar_D_O B_bar_D_O C_bar_D_O D_bar_D_O
global A_A_D B_A_D C_A_D D_A_D
global sigmaProdD_dot_fun
global rSD_goal largeP
global kDOr kDDr kDOv1 alphaDOv kDDv alphaDDv 
global Rjk00 Rjj0

XD=XD(:,indDef);
ND=length(indDef);
Rji0=0.6*R_bar_AD*ones(1,ND);

indAtt=pairDA(indDef);
%rACon=XA(1:2,indAtt);
%vACon=XA(3:4,indAtt);

rDcm=sum(XD(1:2,:),2)/ND;
vDcm=sum(XD(3:4,:),2)/ND;

rA=XA(1:2,:);
vA=XA(3:4,:);
NA=size(XA,2);
NO=size(rVO,2);
for j=1:ND
    F_D=0;  %for adding the fields together
    F_D_dot=[0;0];
    sigmaProd=1; %for product
    sigmaBarProd=1;
    sigmaSum=0;
    sigmaBarSum=0;
    Sigma=zeros(NO+ND-1,1);
    Sigma_dot=zeros(NO+ND-1,1);
    rD_des=XD_des(1:2,j);
    rD=XD(1:2,j);
    vD=XD(3:4,j);
    vD_des=XD_des(3:4,j);
    e_rD=rD-rD_des;
    e_vD=vD-vD_des;
    norm_e_rD=norm(e_rD);
    norm_e_vD=norm(e_vD);
    % check for nearby obstacles
    if(1)
        for k=1:NO
            n=nO(k);
            a=aO(k);
            b=bO(k);
            R_OD=norm(rD-rCO2(1:2,k));
            EmDO=E_m_DO(j,k);
            EbarDO=E_bar_DO(j,k);
            EuDO=E_u_DO(j,k);
            EvDO=E_v_DO(j,k);
            E_DOk=(abs((rD(1)-rCO2(1,k))/a))^(2*n)+(abs((rD(2)-rCO2(2,k))/b))^(2*n)-1;
            E_DOk_dot=2*n*[(abs((rD(1)-rCO2(1,k))/a))^(2*n)/(rD(1)-rCO2(1,k)), (abs((rD(2)-rCO2(2,k))/b))^(2*n)/(rD(2)-rCO2(2,k))]*vD;
            %Modify the boundary elliptic distances
            E_Ddes_Ok=(abs((rD_des(1)-rCO2(1,k))/a))^(2*n)+(abs((rD_des(2)-rCO2(2,k))/b))^(2*n)-1;
            if E_Ddes_Ok<EuDO
                EuDO=E_Ddes_Ok;
                EbarDO=(E_Ddes_Ok-EmDO)/2;
            end
            arr_EDO(j,k)=E_DOk;
            if E_DOk<EuDO
                if E_DOk>EmDO &&  E_DOk<EbarDO
                    sigma=1;
                    sigma_dot=0;
                    sigma_bar=1;
                     FlagDefInObs(j,k)=1;
                elseif  E_DOk>EbarDO &&  E_DOk<EuDO
                    sigma=A_D_O(j,k)*E_DOk^3+B_D_O(j,k)*E_DOk^2+C_D_O(j,k)*E_DOk+D_D_O(j,k);
                    sigma_dot=(3*A_D_O(j,k)*E_DOk^2+2*B_D_O(j,k)*E_DOk+C_D_O(j,k))*E_DOk_dot;
                    sigma_bar=1;
                    if FlagDefInObs(j,k)~=1
                        BetaDv0(j,k)=atan2(vD(2),vD(1)); %Orientation while entering the obstacle field
                       % SpeedD0(j,k)=norm(vD); %Speed while entering the obstacle field
                        mDv0(j,k)=tan(BetaDv0(j,k));
                        beta=atan(sign(-b^(2*n)/(mDv0(j,k)*a^(2*n)))*abs(-b^(2*n)/(mDv0(j,k)*a^(2*n)))^(1/(2*n-1)));
                        rB=((EmDO+1)/(abs(cos(beta))^(2*n)/a^(2*n)+abs(sin(beta))^(2*n)/b^(2*n)))^(1/(2*n));
                        rP0=rCO2(1:2,k)+[rB*cos(beta);rB*sin(beta)];
                        crossProd=cross([cos(BetaDv0(j,k));sin(BetaDv0(j,k));0],[rCO2(1:2,k)-rD;0]);
                        if crossProd(3)>0
                            crossProd2=cross([rP0-rCO2(1:2,k);0],[rCO2(1:2,k)-rD;0]);
                            if crossProd2(3)<0
                                rP0=rCO2(1:2,k)-[rB*cos(beta);rB*sin(beta)];
                            end
                        else
                            crossProd2=cross([rP0-rCO2(1:2,k);0],[rCO2(1:2,k)-rD;0]);
                            if crossProd2(3)>0
                                rP0=rCO2(1:2,k)-[rB*cos(beta);rB*sin(beta)];
                            end
                        end
                        cDv0(j,k)=rP0(2)-mDv0(j,k)*rP0(1);
                    end
                    FlagDefInObs(j,k)=1;
                end
                %Point along the radial direction as the initial condition
                beta=atan2(rD(2)-rCO2(2,k),rD(1)-rCO2(1,k));
                RProj=((EmDO+1)/(abs(cos(beta))^(2*n)/a^(2*n)+abs(sin(beta))^(2*n)/b^(2*n)))^(1/(2*n));
                rP0=[RProj*cos(beta);RProj*sin(beta)]+rCO2(1:2,k);
                f=@(r) [abs(r(1)-rCO2(1,k))^(2*n)/a^(2*n)+abs(r(2)-rCO2(2,k))^(2*n)/b^(2*n)-EmDO-1;...
                        a^(2*n)*sign(r(2)-rCO2(2,k))*abs(r(2)-rCO2(2,k))^(2*n-1)*(rD(1)-r(1))-b^(2*n)*sign(r(1)-rCO2(1,k))*abs(r(1)-rCO2(1,k))^(2*n-1)*(rD(2)-r(2))];
                rDProj(:,k)=fsolve(f,rP0,options);
                %Tangent at the projection point
                beta_barP=atan2(b^(2*n)*sign(rDProj(1,k)-rCO2(1,k))*abs(rDProj(1,k)-rCO2(1,k))^(2*n-1),-a^(2*n)*sign(rDProj(2,k)-rCO2(2,k))*abs(rDProj(2,k)-rCO2(2,k))^(2*n-1));
                rTP=[cos(beta_barP);sin(beta_barP)];
                if rTP'*vD<0  %If the velocity of A is opposite to the tangent then reveres the tangent direction for projection
                    rTP=-rTP;
                end
                %check if the tangent is now crossing the initial velocity
                %direction
                Dv0_hat=[cos(BetaDv0(j,k));sin(BetaDv0(j,k))];
                crossProd=cross([Dv0_hat;0],[cos(beta_barP);sin(beta_barP);0]);
                if crossProd(3)>0
                    rTP=Dv0_hat;
                    rDProj(1,k)=(mDv0(j,k)*rD(2)+rD(1)-mDv0(j,k)*cDv0(j,k))/(1+mDv0(j,k)^2);
                    rDProj(2,k)=mDv0(j,k)*rD(1)+cDv0(j,k);
                end
                vDProj(:,k)=rTP'*vD*rTP;
            elseif E_DOk>EuDO &&  E_DOk<EvDO
                sigma_bar=A_bar_D_O(j,k)*E_DOk^3+B_bar_D_O(j,k)*E_DOk^2+C_bar_D_O(j,k)*E_DOk+D_bar_D_O(j,k);
                sigma=0;
                sigma_dot=0;
                FlagDefInObs(j,k)=0;
            else
                sigma=0;
                sigma_dot=0;
                sigma_bar=0;
                FlagDefInObs(j,k)=0;
            end
            Sigma(k)=sigma;
            Sigma_dot(k)=sigma_dot;
            
            %find vector field direction
%             betaOS=atan2(rD_des(2)-rO(2,k),rD_des(1)-rO(1,k));  %angle between the desired location and the obstacle
%             betaOD=atan2(rD(2)-rO(2,k),rD(1)-rO(1,k));   %angle between the location and the obstacle
%             betaOD_dot=[-(rD(2)-rO(2,k)), (rD(1)-rO(1,k))]*vD/R_OD^2;
%             beta_barOS=atan2(b^(2*n)*sign(rD_des(1)-rO(1,k))*(abs(rD_des(1)-rO(1,k)))^(2*n-1),-a^(2*n)*sign(rD_des(2)-rO(2,k))*(abs(rD_des(2)-rO(2,k)))^(2*n-1));  %angle of the tangent at the safe area
%             
%             beta_bar=atan2(b^(2*n)*sign(rD(1)-rO(1,k))*(abs(rD(1)-rO(1,k)))^(2*n-1),-a^(2*n)*sign(rD(2)-rO(2,k))*(abs(rD(2)-rO(2,k)))^(2*n-1));  %angle of the tangent at the given location
%             if beta_bar<0   % get betaOS between [0,2pi]
%                 beta_bar=beta_bar+2*pi;
%             end
%             Delta_beta=betaOD-betaOS;
%             if Delta_beta<0
%                 Delta_beta=Delta_beta+2*pi;
%             end
%             
%             Delta_beta_bar=beta_barOS-betaOS;
%             if Delta_beta_bar<0
%                 Delta_beta_bar=Delta_beta_bar+2*pi;
%             end
%             if Delta_beta<pi
%                 beta_bar1=beta_bar-Delta_beta_bar+Delta_beta/pi*(Delta_beta_bar-pi);
%                 beta_bar1_dot=(b^(2*n)*(2*n-1)*cos(betaOD)*abs(cos(betaOD))^(2*n-2))/(a^(2*n)*sin(betaOD)*abs(sin(betaOD))^(2*n) + (Delta_beta_bar-pi)/pi)*betaOD_dot;
%             else
%                 beta_bar1=beta_bar-Delta_beta_bar/pi*(Delta_beta-pi);
%                 beta_bar1_dot=(b^(2*n)*(2*n-1)*cos(betaOD)*abs(cos(betaOD))^(2*n-2))/(a^(2*n)*sin(betaOD)*abs(sin(betaOD))^(2*n)-(Delta_beta_bar)/pi)*betaOD_dot;
%             end
%             F_DOk=  [cos(beta_bar1);sin(beta_bar1)];
%             F_D = F_D + sigma*F_DOk;
%             F_D_dot=F_D_dot+sigma_dot*F_DOk+sigma*beta_bar1_dot*[-sin(beta_bar1);cos(beta_bar1)];
            sigmaProd=sigmaProd*(1-sigma);
            sigmaSum=sigmaSum+sigma;
            sigmaBarProd=sigmaBarProd*(1-sigma_bar);
            sigmaBarSum=sigmaBarSum+sigma_bar;
        end
    end
    
    %Potential for projected beta-agents on the static obstacles  
    uDOv=zeros(2,1);
    uDOr=uDOv;
    for k=1:NO
        if Sigma(k)~=0
            Rjk0=Rjk00(j);
            Rj_jk=norm(rD-rDProj(:,k));
            if abs(Rj_jk-R_m_DD)>1e-10
            nabla_rj_VjP= kDOr*Sigma(k)*(rD-rDProj(:,k))/Rj_jk/abs(Rj_jk-R_m_DD)*((Rj_jk-R_m_DD)^2-Rjk0^2)/((Rj_jk-R_m_DD)^2+Rjk0^2);
            else
              nabla_rj_VjP=kDOr*Sigma(k)*(rD-rDProj(:,k))*largeP;
            end
            dv=kDOv1*(vD-vDProj(:,k))*norm(vD-vDProj(:,k))^(alphaDOv-1);
            uDOv=uDOv-Sigma(k)*dv;
            uDOr=uDOr-Sigma(k)*nabla_rj_VjP;
        end
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
                    if R_DjDl>R_m_DD && R_DjDl<R_bar_DD
                        sigma=1;
                        sigma_dot=0;
                    elseif R_DjDl>R_bar_DD && R_DjDl<R_u_DD
                        sigma=A_D_D*R_DjDl^3+B_D_D*R_DjDl^2+C_D_D*R_DjDl+D_D_D;
                        sigma_dot=(3*A_D_D*R_DjDl^2+B_D_D*R_DjDl+C_D_D)*R_DjDl_dot;
                    end
                    Rjj=Rjj0(j);
                    if abs(R_DjDl-R_m_DD)>1e-10
                    nabla_rj_Vjj=kDDr*(-rDjDl)/R_DjDl/abs(R_DjDl-R_m_DD)*((R_DjDl-R_m_DD)^2-Rjj^2)/((R_DjDl-R_m_DD)^2+Rjj^2);
                    else
                        nabla_rj_Vjj=kDDr*(-rDjDl)*largeP;
                    end
              
                    dv=kDDv*(vD-XD(3:4,l))*norm(vD-XD(3:4,l))^(alphaDDv-1);
                    uDDv=uDDv-sigma*dv;
                    uDDr=uDDr-sigma*nabla_rj_Vjj;                    
                else
                    sigma=0;
                    sigma_dot=0;
                end
                Sigma(NO+count)=sigma;
                Sigma_dot(NO+count)=sigma_dot;
%                 F_DjDl=(-rDjDl)/R_DjDl;
%                 F_D = F_D+ sigma*F_DjDl;
%                 F_D_dot=F_D_dot+ sigma_dot*F_DjDl+sigma*(vD-R_DjDl_dot*F_DjDl)/R_DjDl;
                sigmaProd=sigmaProd*(1-sigma);
                sigmaSum=sigmaSum+sigma;
            end
        end
        %sigmaSum
    end
    
    
    %Potential based formation controller
    %Vij(x)=log(k/(x-x0)+(x-x0)/k)
   
%     uDDv=zeros(2,1);
%     uDDr=uDDv;
%     for ii=find(WD(i,:)==1)
%         Rii0=Rii00(i);
%         Ri_ii=norm(rA-XA(1:2,ii));
%         nabla_ri_Vii=kAFr*WA(i,ii)*(rA-XA(1:2,ii))/Ri_ii/abs(Ri_ii-R_m_AA)*((Ri_ii-R_m_AA)^2-Rii0^2)/((Ri_ii-R_m_AA)^2+Rii0^2);
%         dv=kAFv*(vA-XA(3:4,ii))*norm(vA-XA(3:4,ii))^(alphaAFv-1);
%         uDFv=uDFv-dv;
%         uDFr=uDFr-nabla_ri_Vii;
%     end
    
    %Check for nearby attckers to avoid
    ka=1;
    nabla_rj_Vji=zeros(2,1);
    uDA=zeros(2,1);
    for i=1:NA
        rDAi=rA(1:2,i)-rD;
        R_DAi=norm(rDAi);
        R_DAi_dot=(rDAi)'*(vA(1:2,i)-vD)/R_DAi;
        if R_DAi<=rho_c_A
            if R_DAi<R_u_AD
                if R_DAi>R_m_AD && R_DAi<R_bar_AD
                    sigma=1;     sigma_dot=0;
                elseif R_DAi>R_bar_AD && R_DAi<R_u_AD
                    sigma=A_A_D*R_DAi^3+B_A_D*R_DAi^2+C_A_D*R_DAi+D_A_D;
                    sigma_dot=(3*A_A_D*R_DAi^2+2*B_A_D*R_DAi+C_A_D)*R_DAi_dot;
                end
                Rji=0.6*R_bar_AD;
                if abs(R_DAi-R_m_AD)>1e-10
                    nabla_rj_Vji=nabla_rj_Vji+ka*(rD-rA(1:2,i))/R_DAi/abs(R_DAi-R_m_AD)*((R_DAi-R_m_AD)^2-Rji^2)/((R_DAi-R_m_AD)^2+Rji^2);
                else
                    nabla_rj_Vji=nabla_rj_Vji+ka*(rD-rA(1:2,i))*largeP;
                end
                uDA=uDA-sigma*nabla_rj_Vji;
            else
                sigma=0;   sigma_dot=0;
            end
            %             Sigma(NO+j)=sigma;
            %             Sigma_dot(NO+j)=sigma_dot;
            %             F_DAi=(-rDAi)/R_DAi;
            %             F_AD(:,i) = F_AD(:,i) + sigma*F_DAi;
            %             F_AD_dot(:,i) = F_AD_dot(:,i) + sigma_dot*F_DAi+sigma*(vA-R_DAi_dot*F_DAi)/R_DAi;
            %             sigmaProdD=sigmaProdD*(1-sigma);
            %             sigmaSumD=sigmaSumD+sigma;
        end
    end
    
    %for attraction toward the desired locations
%     if norm_e_rD~=0
%         F_DDdes=(-e_rD)/norm_e_rD;
%         norm_e_rD_dot=e_rD'*e_vD/norm_e_rD;
%         F_D = F_D + sigmaProd*F_DDdes;%/Rc^2;
%         s1=num2cell([Sigma',Sigma_dot']);
%         sigmaProd_dot=sigmaProdD_dot_fun(s1{:});
%         F_D_dot=F_D_dot+sigmaProd_dot*F_DDdes+sigmaProd*(-e_vD-norm_e_rD_dot*F_DDdes)/norm_e_rD;
%     end
%     norm_F_D=norm(F_D);
    if (1)% norm_F_D~=0
        %         e_vD_des=-kr1*e_rD*norm_e_rD^(kr2-1);
        %         e_vD_des_dot=kr1*(-norm_e_rD^(kr2-1)*e_vD-(kr2-1)*norm_e_rD^(kr2-3)*e_rD*(e_rD'*e_vD));
%         D_hat=F_D/norm_F_D;
%         norm_F_D_dot=F_D'*F_D_dot/norm_F_D;
%         
%         e_vD_des=kr1*D_hat*norm_e_rD^(kr2);
%         infNorm_e_vD_des=max(abs(e_vD_des(1)),abs(e_vD_des(2)));
%         if infNorm_e_vD_des> v_maxD(j)
%             e_vD_des=e_vD_des*v_maxD(j)/infNorm_e_vD_des;
%         end
%         
%         e_vD_des_dot=kr1*(norm_e_rD^(kr2)*(F_D_dot-D_hat*norm_F_D_dot)/norm_F_D+D_hat*kr2*norm_e_rD^(kr2-2)*(e_rD'*e_vD));
%         
%         duD=C_d*e_vD_des+e_vD_des_dot-kv1*(e_vD-e_vD_des)*(norm(e_vD-e_vD_des))^(kv2-1);
%         
%         duD=C_d*vD_des+e_vD_des_dot-kv1*(vD-e_vD_des)*(norm(vD-e_vD_des))^(kv2-1);
%         
%         e_vD_des=vD-(sigmaBarProd*(vD_des)+sigmaBarSum*kr1*D_hat*norm_e_rD^(kr2));
%         e_vD_des=vD-(sigmaBarProd*(vD_des)+sigmaBarSum*kr1*D_hat*5);
%         
        %         if (norm_e_rD > R_u_AD || abs(-e_rD'*refTraj(1:2,j))<0.7*norm_e_rD*norm(refTraj(1:2,j)))
        %             duD=-0.0005*e_rD*norm_e_rD^(kr2-1)-0.1*vD;
        %         else
%         if (1)%(norm_e_rD > R_u_AD || abs(-e_rD'*refTraj(1:2,j))<0.7*norm_e_rD*norm(refTraj(1:2,j)))
%             % duD=-0.9*e_vD_des*norm( e_vD_des)^(kv2-1)+0.0005*sigmaBarProd*D_hat*norm_e_rD^(kr2)+C_d*(vD-e_vD_des);
%             duD=-0.9*e_vD_des*norm( e_vD_des)^(kv2-1)-0.2*e_rD*norm_e_rD^(kr2-1)-0.1*vD;
%             % duD=-.7*e_vD_des*norm( e_vD_des)^(kv2-1)+uDA;
%             
%         else
%             e_vD_des=vD-refTraj(1:2,j)*norm(vACon(1:2,j))/norm(refTraj(1:2,j));
%             Rji=Rji0(j);
%             Rj_i=norm(rACon(1:2,j)-rD);
%             R_DAi=Rj_i;
%             uDA= ka*(rD-rACon(1:2,j))/R_DAi/abs(R_DAi-R_m_AD)*((R_DAi-R_m_AD)^2-Rji^2)/((R_DAi-R_m_AD)^2+Rji^2); %just for checking, original control is defined above
%             %duD = -.7*e_vD_des*norm( e_vD_des)^(kv2-1)+uDA;
%             duD = -.005*(rD-rS)+uDA;
%             
%         end
        %duD = -0.1*(rD-rSD_goal(:,j))*tanh(norm(rD-rSD_goal(:,j)))/norm(rD-rSD_goal(:,j))+uDOv+uDOr+uDDv+uDDr;
         %duD = -0.2*(rD-rD_des)*tanh(norm(rD-rD_des))/norm(rD-rD_des)+uDOv+uDOr+uDDv+uDDr;
         duD = -0.2*(rD-rD_des)*norm(rD-rD_des)^(0.5-1)+uDOv+uDOr+uDDv+uDDr;
        %duD=C_d*vD_des-kv1*(vD-e_vD_des)*(norm(vD-e_vD_des))^(kv2-1);
        uD0=duD;%+C_d*vD_des;
        uD_des0=uD0;%norm(uD0)*F_D/norm_F_D;
    else
        uD_des0=ones(2,1);
    end
    
    uD(:,j)=uD_des0;
    
    
    %norm_uD=norm(uD(:,j));
    infNorm_uD=max(abs(uD(1,j)),abs(uD(2,j)));
    if infNorm_uD>u_maxD(j)
        uD(:,j)=uD(:,j)*u_maxD(j)/infNorm_uD;
    end
end
end
