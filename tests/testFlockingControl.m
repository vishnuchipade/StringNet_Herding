

AllParameters;
% ------------------------------
% Author: Vishnu S. Chipade
% Email: vishnuc@umich.edu
% ------------------------------

%ND=0;
t=0;
Time(1)=t;
t_max=1000; endGame=0;
Niter=floor(t_max/dt);
%Allocate memory
X=zeros(4*(ND+NA),Niter+1);
for i=1:NA
    vA_dot(:,i)=[0;0];
end
XA=XA0;
XD=XD0;
X(:,1)=[XA(:);XD(:)];
SD(:,1)=zeros(ND,1);
SD_arr(:,1)=SD;

%Communication topology of the attackers and formation distances
global Rii_bar
Rii0=RA0*sqrt(2*(1-cos(2*pi/NA)));
Rii1=Rii0*sqrt(2*(1-cos(180-2*pi/NA)));
Rii_bar=Rii1*ones(NA);
if NA>1
    WA=ones(NA);
    NNA=2;
    %NAC=NA-1;  $One center agent, others on the boundary
    NAC=NA;  % all at the periphery, No center agent
    temp=[1:NAC,1:NAC];
    for i=1:NAC
        WA(i,temp([i+1:i+NNA/2,i+NAC-NNA+1:i+NAC-NNA/2]))=ones(1,NNA);
        Rii_bar(i,temp([i+1:i+NNA/2,i+NAC-NNA+1:i+NAC-NNA/2]))=Rii0;
        Rii_bar(i,i)=0;
        WA(i,i)=0;
    end
    %WA(NA,:)=1;  WA(:,NA)=1; WA(NA,NA)=0;
else
    WA=0;
end

[WA, Rii_bar]=findCommGraphAndFormDist(NA,0,RA0);
rCO2=obs.rCO2;

dpsiAi0=0.02;
vA_ref_max=4;
global dpsi_var
syms dpsi_var
flagAttInObs=zeros(NA,NO);
FlagDefInObs=zeros(ND,NO);
betaAv0=zeros(NA,NO);
betaDv0=zeros(ND,NO);
mDv0=betaDv0;
cDv0=mDv0;
mAv0=betaAv0;
cAv0=mAv0;
speedA0=betaAv0;
assignment=[1:ND];
distTol=3;
countDDes=0;
for ti=1:Niter
    if mod(ti,1000)==0
        ti
    end
    %     if ti==1728
    %         ti
    %     end
    %desired positions for the leader attackers
    for i=1:NA
        thetaA=360-360*(i-1)/NA;
        XA_goal(1:2,i)=rP+[RA0*cos(thetaA*pi/180), RA0*sin(thetaA*pi/180);]';
        
        %For flocking
        XA_goal(1:2,i)=rP;
        XA_goal(3:4,i)=[0,0]';
    end
    %rD=rD_des;
    
    %[uA, uA0, vA_des, vA_des_dot,rAProj1,flagAttInObs,betaAv0,speedA0,mAv0,cAv0, sigmaProdD, FA, FA_dot] = controlAttacker3(XA,XA_goal,vA,flagAttInObs,betaAv0,speedA0,mAv0,cAv0,XD,WA,WDString,NA,ND);
    %    vA_Des(:,ti)=vA_des(:);
    %     vA=[0,0]';
    %     vA_dot=[0,0]';
    
    %desired locations of the defenders for herding
    %[Psi, Psi_dot, Psi_ddot] = formationRefTraj(XA(1:2,:),XA(3:4,:),uA(:,:),rO,rS);
    % [psi, psi_dot, psi_ddot] = formationRefTraj(rA,vA,vA_dot,rO,rS);
    %Psi_des0(ti)=Psi;
    Psi=zeros(NA,1);
    Psi_dot=Psi;
    Psi_ddot=Psi;
    
    %check if the formation has entered the safe area
    countAInS=0;
    for ii=1:NA
        if norm(rA(:,ii)-rS)<rho_S
            countAInS=countAInS+1;
        end
    end   
 
    
    %%
    [uA, uA0, vA_des, vA_des_dot,rAProj1,flagAttInObs,betaAv0, speedA0,mAv0,cAv0, SigmaProdD, FA, FA_dot] = controlAttacker3(XA,XA_goal,vA,flagAttInObs,betaAv0,speedA0,mAv0,cAv0,[],WA,[],NA,0);
    rAProj_arr(:,:,ti)=rAProj1;  %row=attacker, column=obstacle, slice=time
    SigmaProdD_arr(:,ti)=SigmaProdD;
    uA(:,NA)=[0,0]';
    %Control action for the defenders
    %[uD, arr_EDO,e_vD_des]= controlDefender4(XD,indDef, XD_des, XD_des_dot, XA,pairDA,refTraj, ND);
    rAcm=sum(XA(1:2,:),2)/NA;
    vAcm=sum(XA(3:4,:),2)/NA;
    rAcm_arr(:,ti)=rAcm;
    vAcm_arr(:,ti)=vAcm;
    
    vA_dot_arr(:,ti)=vA_dot(:);
    FA_dot_arr(:,ti)=FA_dot(:);
    FA_arr(:,ti)=FA(:);
    if norm(vA_dot)>10
        norm(vA_dot);
    end
    %Psi(ti)=atan2(vA(2),vA(1));
    uA(:,NA)=[0,0]';
    U(:,ti)=[uA(:);zeros(2*ND,1)];
  
    X(:,ti+1)=modifiedDIDynamics(X(:,ti),U(:,ti));%,NA,ND);
    XA=reshape(X(1:4*NA,ti+1),[4,NA]);
    rA=XA(1:2,:);
    vA=XA(3:4,:);
  
 
    
    %Find critical distances for the attackers
    for i=1:NA-1
        for l=i+1:NA
            RAiAl(l)=norm(rA(:,i)-rA(:,l));
        end
        arr_minRAA(i)=min(RAiAl(i+1:end));
        %super-elliptic distance of the attacker from the obstacles
        for kk=1:NO
            n=nO(kk);
            a=aO(kk);
            b=bO(kk);
            % RDjOl(l)=norm(rD(:,j)-rO(1:2,l));
            EAOl(kk)=(abs((rA(1,i)-rCO2(1,kk))/a))^(2*n)+(abs((rA(2,i)-rCO2(2,kk))/b))^(2*n)-1;
            EAO_rel(i,kk)=E_m_AO(i,kk)/EAOl(kk);
        end
    end
    if ~isempty(arr_minRAA)
        minRAA(ti+1)=min(arr_minRAA);
    end
    if NA>1
        minEAO(ti+1)=max(max(EAO_rel));
    end
    
    
    t=t+dt;
    time(ti+1)=t;
    %ti=ti+1;
end 