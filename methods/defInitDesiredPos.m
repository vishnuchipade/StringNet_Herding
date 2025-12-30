function [defDesForm, motionPlan]=defInitDesiredPos(XA0,XD,NA,ND,tanG,RD0,v_maxA,DeltaT_s)
% ------------------------------
% Author: Vishnu S. Chipade
% Email: vishnuc@umich.edu
% ------------------------------

global rP rho_sn rho_P v rho_Acon rho_safe
%Initial reference semicircular formation
rAcm=sum(XA0(1:2,:),2)/NA;
vAcm=sum(XA0(3:4,:),2)/NA;
%thetaAcm0=atan2(vAcm(2),vAcm(1))+20*pi/180;   %Added extra 20 just to check
thetaAcm0=atan2(rAcm(2)-rP(2),rAcm(1)-rP(1))+20*pi/180;   %Added extra 20 just to check
if thetaAcm0<0
    thetaAcm0=thetaAcm0+2*pi;
end
[path_Acm,tanG_prime_Acm,rE_dash] = findShortestPath(tanG,rAcm,rP);
for i=2:length(path_Acm.nodes)-1
    ang(i)=atan2(path_Acm.rV(2,i)-path_Acm.rVC(2,i),path_Acm.rV(1,i)-path_Acm.rVC(1,i));
    if ang(i)<0
        ang=ang+2*pi;
    end
    path_Acm.rV(:,i)=path_Acm.rVC(:,i)+(RD0+rho_safe)*[cos(ang(i)),sin(ang(i))]';
    if mod(i,2)==0
        path_Acm.P(i-1)=norm(path_Acm.rV(:,i)-path_Acm.rV(:,i-1));
    else
         path_Acm.P(i-1)=(RD0+rho_safe)*(-2*path_Acm.segType(i-1)+3)*(ang(i)-ang(i-1));
    end
    path_Acm.S(i)=path_Acm.S(i-1)+path_Acm.P(i-1);
end
path_Acm.P(i)=norm(path_Acm.rV(:,i)-path_Acm.rV(:,i+1));
path_Acm.S(i+1)=path_Acm.S(i)+path_Acm.P(i);

Qa=path_Acm.S(end);
q1=0;
q2=Qa-rho_P;


f=10000;
i=0;
tComp=0;
df=10000;
while abs(f)>3 && abs(df)>1
    i=i+1;
    q=(q1+q2)/2;
    [rDFc0, thetaAcm0]=findCoordOnPath(q,path_Acm);
    
    %rDFc0=rP+(3*rho_P+rho_sn)*[cos(thetaAcm0), sin(thetaAcm0)]';
    thetaD=thetaAcm0-pi/2;
    for j=1:ND
        %RD0=rho_sn;%3*rDmin;
        
        rD_des(:,j)=rDFc0+0.5*(ND-1)*RD0*[cos(thetaD), sin(thetaD);]'-RD0*[cos(thetaD), sin(thetaD);]'*(j-1);
        XD_des0(:,j)=[rD_des(:,j);0;0];
        XD_des_dot0(:,j)=[0,0,0,0]';
    end
    
    motionPlan=motionPlanForDefOpenForm(tanG,XD,XD_des0,ND,0,1);
    tComp=tComp+motionPlan.tComp;
    df=f-(motionPlan.maxOptT-(q)/v_maxA+DeltaT_s);
    f=motionPlan.maxOptT-(q)/v_maxA+DeltaT_s;
    if f<0
        q2=q;
    else
        q1=q;
    end
    F(i)=f;
    Q(i)=q;
end
motionPlan.tComp=tComp;
defDesForm.rDFc0=rDFc0;
defDesForm.XD_des0=XD_des0;
defDesForm.XD_des_dot0=XD_des_dot0;
defDesForm.phi=thetaAcm0;
defDesForm.phi_dot=0;
end