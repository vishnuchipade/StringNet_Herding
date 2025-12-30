function [XD_des, XD_des_dot,dpsiAi0] = desiredPosDefender(XD,indDef,XA,pairDA,uA,uA0,refTraj,Psi,Psi_dot,Deltaj,dpsiAi0,NA,ND)
% ------------------------------
% Author: Vishnu S. Chipade
% Email: vishnuc@umich.edu
% ------------------------------

%This function calculates the desired position of the given defenders around the
%attackers based on the mapping of defenders with the attckers and the
%desired reference direction
global options dpsi_var R_m_AD R_bar_AD R_u_AD Rij0
global kADr
rD=XD(1:2,indDef);
vD=XD(3:4,indDef);

%get the states of all the attackers and find their CoM
rA=XA(1:2,:);

vA=XA(3:4,:);
rAcm=sum(rA,2)/NA;
vAcm=sum(vA,2)/NA;
vAcm_dot=sum(uA,2)/NA;
for i=1:NA
dist(i)=norm(rA(:,i)-rAcm);
end
max_dist=max(dist);
%get the positions and the velocities of the attackers to be controlled in
%the order such that first attacker in the array corresponds to the first
%defender as per the list indDef
indAtt=pairDA(indDef);
rAcon=XA(1:2,indAtt);
vAcon=XA(3:4,indAtt);
Psi=Psi(indAtt);
Psi_dot=Psi_dot(indAtt);
%Desired trajectory of the CoM of the attackers
% vAc_des=refTraj(1:2,1);
% psiAc_des=atan2(vAc_des(2),vAc_des(1));
% vAc_des_dot=refTraj(1:2,2);
%psi_ddot=refTraj(3,1);\
vA_des=refTraj(1:2,:);%
vA_des_dot=refTraj(3:4,:);

% uAc=sum(uA,2)/NA;
% uAc0=sum(uA0,2)/NA;
% duAc=-uAc0-(vAc-vAc_des)+vAc_des_dot;
% norm_duAc=norm(duAc);
% psi_duAc=atan2(duAc(2),duAc(1));
% RAD=ND/norm_duAc;
% 
% mat=zeros(2,2);
% for i=1:length(indDef)
%     psi_rA_rAc(i)=atan2(rAc(2)-rAcon(2,i),rAc(1)-rAcon(1,i));
%     mat=mat+[cos(psi_rA_rAc(i)), -sin(psi_rA_rAc(i));sin(psi_rA_rAc(i)) , cos(psi_rA_rAc(i))];
% end
% 
%  var=(norm_duAc*mat\duAc)/norm(norm_duAc*mat\duAc);
%  ang=atan2(var(2),var(1));
% f=@(dpsi_var) tan(psi_duAc)-tan((mat(2,1)*sin(dpsi_var)+mat(2,2)*cos(dpsi_var))/(mat(1,1)*cos(dpsi_var)+mat(1,2)*sin(dpsi_var)));
% 
% dpsiAi=fsolve(f,dpsiAi0,options);
%dpsiAi=ang;
%dpsiAi0=dpsiAi;
for j=1:length(indDef)
    ji=j;
    rAi=rA(:,ji);
    vAi=vA(:,ji);
    vAi_dot=uA(:,ji);
    vAi_des=vA_des(1:2,ji);
    vAi_des_dot=vA_des_dot(1:2,ji);
    if norm(rAi-rD(:,ji))> 1*R_u_AD
    duAi=-uA0(:,ji);%-(vAi-vAi_des)+vAi_des_dot;
    else
    duAi=-uA0(:,ji)-(vAi-vAi_des)+vAi_des_dot;
    end
   norm_duAi=norm(duAi);
   
   % RAD=max_dist+(R_u_AD+norm_duAi*R_m_AD)/(norm_duAi+1);
    %RAD=(R_u_AD+norm_duAi*R_m_AD)/(norm_duAi+1);
   % RAD=max_dist+R_bar_AD;
   %RAD=0.6*R_bar_AD;
    
   if ~isnan(norm_duAi)
    RAD0=roots([norm_duAi,norm_duAi*Rij0(ji)^2,kADr,-kADr*Rij0(ji)^2]);
   else
       RAD0=0.1;
   end
    RAD0=RAD0(imag(RAD0)==0);
    if max(RAD0)>0
    RAD=max(RAD0)+R_m_AD;
    else
        disp('Negative radius');
        break;
    end
    
    
    %psiAi_des=Psi(i);
    psiAi_des=atan2(duAi(2),duAi(1));
   % rD_des(:,j)=rAi+RAD*[cos(psi_rA_rAc(j)+dpsiAi+pi), sin(psi_rA_rAc(j)+dpsiAi+pi)]';  %constant radius from attacker
 %  rD_des(:,j)=rAcm+RAD*[cos(psiAi_des+pi+Deltaj(j)), sin(psiAi_des+pi+Deltaj(j))]';
%     rD_des(:,j)=rAcm+RAD*[cos(psiAi_des+pi), sin(psiAi_des+pi)]';
  %  vD_des(:,j)=vAcm;%+RAD*[-sin(psi_prime+pi+Deltaj(j)), cos(psi_prime+pi+Deltaj(j))]'*psi_prime_dot;
  %   aD_des(:,j)=vAcm_dot;%+RAD*[-sin(psi_prime+pi+Deltaj(j)), cos(psi_prime+pi+Deltaj(j))]'*psi_prime_ddot+...
    %         +RAD*[cos(psi_prime+pi+Deltaj(j)), sin(psi_prime+pi+Deltaj(j))]'*psi_prime_dot^2;
    
    %rD_des(:,j)=rAi+RAD*[cos(psiAi_des+pi+Deltaj(j)), sin(psiAi_des+pi+Deltaj(j))]';
    
     rD_des(:,j)=rAi+RAD*[cos(psiAi_des+pi), sin(psiAi_des+pi)]';
     vD_des(:,j)=vAi;%+RAD*[-sin(psi_prime+pi+Deltaj(j)), cos(psi_prime+pi+Deltaj(j))]'*psi_prime_dot;
     aD_des(:,j)=vAi_dot;
    
    XD_des0(:,j)=[rD_des(:,j);vD_des(:,j)];
    XD_des_dot0(:,j)=[vD_des(:,j);aD_des(:,j)];
end
XD_des(:,indDef)=XD_des0;
XD_des_dot(:,indDef)=XD_des_dot0;
end

