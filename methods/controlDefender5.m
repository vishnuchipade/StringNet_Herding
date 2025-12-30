function [ uD, arr_EDO] = controlDefender5(XD, SD, indDef, assign, XD_des, XD_des_dot, Path, pathVel, t, startTime, ND )
% ------------------------------
% Author: Vishnu S. Chipade
% Email: vishnuc@umich.edu
% ------------------------------

%Time optimal control along the shortest path
global rho_safe rS obs nO aO bO kr0 kr1 kr2 kv1 kv2 RD_con C_d rho_c_A options
global E_m_DO E_bar_DO E_u_DO E_v_DO R_m_DD R_bar_DD R_u_DD R_u_AD R_m_AD R_bar_AD
global A_D_O B_D_O C_D_O D_D_O A_D_D B_D_D C_D_D D_D_D
global A_bar_D_O B_bar_D_O C_bar_D_O D_bar_D_O
global A_A_D B_A_D C_A_D D_A_D
global sigmaProdD_dot_fun
global rSD_goal largeP
global kDOr kDDr kDOv1 alphaDOv kDDv alphaDDv
global Rjk00 Rjj0

global v_maxDC v_maxD u_maxD

XD=XD(:,indDef);
XD_des=XD_des(:,assign);
%ND=length(indDef);
Rji0=0.6*R_bar_AD*ones(1,ND);

%indAtt=pairDA(indDef);
%rACon=XA(1:2,indAtt);
%vACon=XA(3:4,indAtt);

rDcm=sum(XD(1:2,1:ND),2)/ND;
vDcm=sum(XD(3:4,1:ND),2)/ND;

rD=XD(1:2,:);
vD=XD(3:4,:);
% NA=size(XA,2);
NO=obs.NO;
rCO2=obs.rCO2;
rVO=obs.rVO;

uD(:,ND+1)=[0,0]';

for j=1:ND
    if t>=startTime(j)
        s_ddot=0;
        s_dot=norm(vD(:,j));
        S=Path{j}.S;
        P=Path{j}.P;
        rV=Path{j}.rV;
        rVC=Path{j}.rVC;
        segType=Path{j}.segType;
        v_bar=pathVel{j}.v_bar;
        s_bar1=pathVel{j}.s_bar1;
        s_bar2=pathVel{j}.s_bar2;
        ind0=find(S<=SD(j));
        indS=ind0(end);
        if indS<=length(P)
            %         if indS==length(P)
            %             1;
            %         end
            %         if indS==1
            %             S_prev=0;
            %         else
            %             S_prev=S(indS-1);
            %         end
            S_prev=S(indS);
            %If straight segment
            if mod(indS,2)==1
                if  SD(j)<s_bar1(indS)+S_prev  %s_dot<v_bar(indS) &&
                    s_ddot=u_maxD(j);%-C_d*s_dot^2;
                elseif  SD(j)>s_bar2(indS)+S_prev  %s_dot<v_bar(indS) &&
                    s_ddot=-(u_maxD(j));%-C_d*s_dot^2;
                else
                    s_ddot=u_maxD(j);
                end
                uD(:,j)=s_ddot*(rV(:,indS+1)-rV(:,indS))/norm(rV(:,indS+1)-rV(:,indS));
            end
            
            %If circular arc segment
            if mod(indS,2)==0
                if  SD(j)<s_bar1(indS)+S_prev  %s_dot<v_bar(indS) &&
                    s_ddot=(u_maxD(j)-s_dot^4/(rho_safe^2*u_maxD(j))); %-C_d*s_dot^2+
                elseif  SD(j)>s_bar2(indS)+S_prev  %s_dot<v_bar(indS) &&
                    s_ddot=-(u_maxD(j)-s_dot^4/(rho_safe^2*u_maxD(j)));  %-C_d*s_dot^2
                else
                    s_ddot=0;
                end
                dr=XD(1:2,j)-rVC(:,indS);
                dr_hat=dr/norm(dr);
                uD(:,j)=s_ddot*(-1)^(segType(indS)-1)*[- dr_hat(2); dr_hat(1)]+s_dot^2/rho_safe*(-dr_hat);
            end
            if norm(uD(:,j))>100
                1
            end
        else
            uD(:,j)=-(rD(1:2,j)-XD_des(1:2,j))-vD(1:2,j);
        end
        
    else %startTime hasn't crossed (stay stationary)
        uD(:,j)=[0,0]';
    end
    %Find the distance to obstacles
    for k=1:NO
        E_DOk=Inf;
        for kk=1:length(rVO{k}(1,:))-1
            dx=rVO{k}(1,kk+1)-rVO{k}(1,kk);
            dy=rVO{k}(2,kk+1)-rVO{k}(2,kk);
            mDD=dy/dx;  %slope of line joining the defenders
            cDD=rVO{k}(2,kk)-mDD*rVO{k}(1,kk);
            rDProjO(1,1)=(mDD*rD(2)+rD(1)-mDD*cDD)/(1+mDD^2);
            rDProjO(2,1)=mDD*rDProjO(1,1)+cDD;
            if mDD<1e16
                lambdaDP=(rDProjO(1,1)-rVO{k}(1,kk))/dx;
            else
                lambdaDP=(rDProjO(2,1)-rVO{k}(2,kk))/dy;
            end
            if lambdaDP<0
                rDProjO(:,1)=rVO{k}(1:2,kk);
            elseif lambdaDP>1
                rDProjO(:,1)=rVO{k}(1:2,kk+1);
            end
            E_DOk=min(E_DOk,norm(rD-rDProjO));
        end
        E_DOk=E_DOk-rho_safe;
        arr_EDO(j,k)=E_DOk;
    end
    
end
