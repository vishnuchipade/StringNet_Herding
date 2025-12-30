function [ uD, arr_EDO] = controlDefenderFormation2(XD, indDef, XD_des, WD, WDString, Rjj_bar, XA,NA,ND)
% ------------------------------
% Author: Vishnu S. Chipade
% Email: vishnuc@umich.edu
% ------------------------------

%Flocking type control for defenders
global obs rS nO aO bO kr0 kr1 kr2 kv1 kv2 RD_con v_maxD u_maxD u_maxA C_d rho_c_A options
global E_m_DO E_bar_DO E_u_DO E_v_DO R_m_DD R_bar_DD R_u_DD R_u_AD R_m_AD R_bar_AD
global A_D_O B_D_O C_D_O D_D_O A_D_D B_D_D C_D_D D_D_D
global A_bar_D_O B_bar_D_O C_bar_D_O D_bar_D_O
global A_A_D B_A_D C_A_D D_A_D
global sigmaProdD_dot_fun
global rSD_goal largeP
global kDOr kDDr kDOv1 kDOv2 alphaDOv kDDv alphaDDv kDFr kDFr2 kDFv alphaDFr alphaDFv kDRr kDRv
global Rjk00 
global rho_Acon rho_sn

global R_m2_DD 

XD=XD(:,indDef);
XD_des=XD_des(:,indDef);

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
XD=[XD(:,1:ND),[rAcm;vAcm]];

uD(:,ND+1)=[0,0]';

for j=1:ND
    rD=XD(1:2,j);
    vD=XD(3:4,j);
    
    %For converging to the formation
    %Vij(x)=log(k/(x-x0)+(x-x0)/k)
    
    uDFv=zeros(2,1);
    uDFr=uDFv;
    for jj=find(WD(j,:)==1)
        Rjj0=Rjj_bar(j,jj);
        Rj_jj=norm(rD-XD(1:2,jj));
        if jj<=ND
            flagString=WDString(j,jj);
        else           
            flagString=1;
        end
        
        if flagString~=1 
        %for Vij(x)=log(k/(x-x0)+(x-x0)/k)        
        nabla_rj_Vjj=kDFr*WD(j,jj)*(rD-XD(1:2,jj))/Rj_jj/abs(Rj_jj-R_m_DD)*((Rj_jj-R_m_DD)^2-Rjj0^2)/((Rj_jj-R_m_DD)^2+Rjj0^2); 
        else
        %for Vij(x)=log((c*x-k)/(x-xm)+(x-xm)/(c*x-k));  (xm=R_m_DD)
        c=(Rjj0-R_m_DD)/(R_m2_DD-Rjj0);
        k=R_m2_DD*c;      
        temp1=(c^2-1)*Rj_jj^2-2*c*k*Rj_jj+k^2-R_m_DD^2+2*R_m_DD*Rj_jj;
        temp2=(c^2+1)*Rj_jj^2-2*c*k*Rj_jj+k^2+R_m_DD^2-2*R_m_DD*Rj_jj;
        nabla_rj_Vjj=kDFr*WD(j,jj)*(rD-XD(1:2,jj))/Rj_jj*(-(k-c*R_m_DD)/((R_m_DD-Rj_jj)*(c*Rj_jj-k))*temp1/temp2);
        end
        %if jj==ND+1 && Rj_jj<rho_sn
            dv=kDFv*(vD-XD(3:4,jj))*norm(vD-XD(3:4,jj))^(alphaDFv-1);
        %else
            %dv=[0,0]';
        %end
        uDFv=uDFv-dv;
        uDFr=uDFr-nabla_rj_Vjj;
    end
        
    %control corresponding to the desired location
    
   
    %final control
    if norm(rD-XD_des(1:2,j))>rho_sn
   uD(:,j)=uDFv+uDFr-kDRr*(rD-XD_des(1:2,j));
    else
        uD(:,j)=uDFv+uDFr-kDRr*(rD-XD_des(1:2,j))-kDRv*(vD-XD_des(3:4,j));
    end
   
   
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

