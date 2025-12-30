%This script simulates herding of group of attackers by a swarm of
%defenders

% ------------------------------
% Author: Vishnu S. Chipade
% Email: vishnuc@umich.edu
% ------------------------------

setupPaths  % <-- this adds all paths

clearvars;% -except tanG Path pathVel

flagExp=1;
flagGazeboExp=0;
if flagExp==1
    AllParametersExperiment;
elseif flagGazeboExp==1
    AllParametersGazeboExperiment;
else
    AllParameters;
end
%AllParametersExperiment;

symDerivative;
%psi_prime_dot_fun=matlabFunction(psi_prime_dot,'Vars',[Del,FAO, pp,p,p_dot,g,g_dot]);
global options
options = optimset('Display','off','MaxIter',1000);


t=0;
Time(1)=t;
t_max=2000; endGame=0;
Niter=floor(t_max/dt);
%Allocate memory
X=zeros(4*(ND+NA+1),Niter+1);
for i=1:NA
    vA_dot(:,i)=[0;0];
end
XA=XA0;
XD=XD0;

%measurements
Cov_YA=diag([0,0,0,0]);
%Cov_YA=diag([.8,.8,.5,0.5]);
YA=XA+mvnrnd(zeros(4,1),Cov_YA,NA)';

SD(:,1)=zeros(ND,1);
SD_arr(:,1)=SD;
rAcm=sum(YA(1:2,:),2)/NA;
rAcm=sum(YA(3:4,:),2)/NA;
rDcm=sum(XD(1:2,1:ND),2)/ND;
vDcm=sum(XD(3:4,1:ND),2)/ND;


flagNDeven=mod(ND,2);

%Initial control
uA=vA;

%Initialization of vectors and matrices
RPA_arr=zeros(1,Niter+1);
RAD_arr=zeros(1,Niter+1);
Psi=zeros(1,Niter+1);
Psi_prime=zeros(1,Niter+1);
Psi_des=zeros(1,Niter+1);
Psi_des_dot=zeros(1,Niter+1);
Psi_des_ddot=zeros(1,Niter+1);
Psi_prime=zeros(1,Niter+1);
Psi_prime_dot=zeros(1,Niter+1);
Psi_prime_ddot=zeros(1,Niter+1);
vA_dot_arr=zeros(2*NA,Niter+1);
FA_dot_arr=zeros(2*NA,Niter+1);
XD_Des=zeros(4*ND,Niter+1);
U=zeros(2*(ND+NA+1),Niter+1);
time=zeros(1,Niter+1);
minEAO=zeros(1,Niter+1);
minEDO=minEAO;
minRDD=minEAO;
minRAD=minEAO;
minRAA=minRAD;
arr_norm_F_AO=minEAO;

%Other initializations
AOut=[];
NAOut=0;
AIcount=0;  %counter for Attacker interceptions;
Winner='None';

safe_flag=0;
RPA=norm(rA-rP);
epsilon=pi/100;
RAD=RAD_max;

RA0_des=0;
%Communication topology of the attackers and formation distances
global Rii_bar
Rii0=RA0*sqrt(2*(1-cos(2*pi/NA)));
Rii1=Rii0*sqrt(2*(1-cos(180-2*pi/NA)));
[WA, Rii_bar]=findCommGraphAndFormDist(NA,1,RA0);  %second argument is shapeId:

%Neighboring graph for the defenders during the closed StringNet formation
[WD_closed, Rjj_bar_closed]=findCommGraphAndFormDist(ND+1,2,rho_sn);

%Neighboring graph for the defenders while moving the open StringNet formation
%around the attackers
if flagExp==1
    RDF_open=3.4;
elseif flagGazeboExp==1
    RDF_open=1.5*rho_sn;
else
    RDF_open=1.5*rho_sn;
end
[WD_open, Rjj_bar_open]=findCommGraphAndFormDist(ND+1,4,RDF_open);

R_DD_string=1.4*Rjj_bar_closed(1,2);

RDF_closed=1.1*rho_sn;

%Initial string chain
% if ND>1
%     WDString=zeros(ND);
%     NND=2;
%     NAC=ND;
%     temp=[1:NAC,1:NAC];
%     %Open chain
%     for j=2:NAC-1
%         WDString(j,temp([j+1:j+NND/2,j+NAC-NND+1:j+NAC-NND/2]))=ones(1,NND);
%     end
%     WDString(1,2)=1;
%     WDString(NAC,NAC-1)=1;
% else
%     WDString=0;
% end
WDString=zeros(ND);
flagHerd=0;
flagDForm=0;

%%
%C1-Tangent graph
flagTanGraph=0;
if(flagTanGraph)
    tanG=tangentGraph(rVO);
    E_m=0;
    for k=1:NO
        rVO2=tanG.rVO2{k};
        rVO_temp=tanG.rVO{k};
        rVO_temp=rVO_temp(:,1:end-1);
        nVO=tanG.nVO(k);
        polyin=polyshape(rVO_temp(1,:),rVO_temp(2,:));
        [rCO2(1,k),rCO2(2,k)]=centroid(polyin);
        %distances between the centroid and the line segments
        for ki=1:2*nVO
            if ki~=2*nVO
                r1=rVO2(:,ki);
                r2=rVO2(:,ki+1);
            else
                r1=rVO2(:,ki);
                r2=rVO2(:,1);
            end
            dx=r2(1)-r1(1);
            dy=r2(2)-r1(2);
            distL(ki)=abs(dy*rCO2(1,k)-dx*rCO2(2,k)+r2(1)*r1(2)-r2(2)*r1(1))/(sqrt(dx^2+dy^2));
        end
        
        %distances between the centroid and the vertices
        for ki=1:nVO
            distV(ki)=norm(rVO{k}(:,ki)-rCO2(:,k));
        end
        
        %E_m0=pi/2*max(distV)/min(distL);
        E_m0=(tanG.PeriO(k)/2)/(2*min(distL));
        %factor for ellipse filter
        if E_m0>E_m
            E_m=E_m0;
        end
    end
    tanG.rCO2=rCO2;
    tanG.E_m=E_m;
    obs.rCO2=rCO2;
    if flagExp==1
        save('tempDataExp.mat','obs','tanG','-append');
    elseif flagGazeboExp==1
        save('tempDataGazeboExp.mat','obs','tanG','-append');
    else
        save('tempData.mat','obs','tanG','-append');
    end
else
    if flagExp==1
        load('tempDataExp.mat');
    elseif flagGazeboExp==1
        load('tempDataGazeboExp.mat')
    else
        load('tempData.mat');
    end
    flagTanGraph=0;
    rCO2=obs.rCO2;
end

%Find the motion plan for defenders to achieve open formation
flagMotionFlag=0;
if flagExp==1
    Delta_Ts=75;
elseif flagGazeboExp==1
    Delta_Ts=16;
else
    Delta_Ts=20;
end
if flagMotionFlag
    %Find the desired initial open formation for the defenders to coverge to
    [defDesForm, motionPlan]=defInitDesiredPos(YA,XD0,NA,ND,tanG,RDF_open,v_maxA(1),Delta_Ts);
    XD_des0=defDesForm.XD_des0;
    %motionPlan=motionPlanForDefOpenForm(tanG,XD0,XD_des0,ND);
    if flagExp==1
        save('tempDataExp.mat','motionPlan','defDesForm','-append');
    elseif flagGazeboExp==1
        save('tempDataGazeboExp.mat','motionPlan','defDesForm','-append');
    else
        save('tempData.mat','motionPlan','defDesForm','-append');
    end
end

rDFc_des0=defDesForm.rDFc0;
XD_des0=defDesForm.XD_des0;
XD_des_dot0=defDesForm.XD_des_dot0;
phi=defDesForm.phi+pi;
phi_dot=defDesForm.phi_dot;
clear defDesForm;

%Extract Motion plan
assign=motionPlan.assign;
assignCost=motionPlan.assignCost;
maxOptT=motionPlan.maxOptT;
assignedPath=motionPlan.Path;
assignedTanG_prime=motionPlan.tanG_prime;
assignedInterSec=motionPlan.interSec;
assignedPathVel=motionPlan.pathVel;
assignedOptT=motionPlan.optT;
leadTime=motionPlan.leadTime;
startTime=motionPlan.startTime+ones(ND,1)*25;
clear motionPlan;
clear tanG;

%XD=XD_des0(:,assign);   %for checking the rotational controller
%add virtual defender
XD=[XD,[rDFc_des0;0;0]];
X(:,1)=[XA(:);XD(:)];

%save('tempData.mat','obs','tanG','motionPlan','defDesForm');

%Defenders indices for the formation
indDef=[1:ND+1]';
na=length(assign);
nid=length(indDef);
indDef([assign;na+1:nid])=indDef;

%%
%minimum distance between the defenders
if ND>1
    for j=1:ND-1
        for l=j+1:ND
            RDjDl(l)=norm(rD(:,j)-rD(:,l));
        end
        arr_minRDD(j)=min(RDjDl(j+1:end));
        arr_minRAD(j)=norm(rD(:,j)-rA);
    end
    minRDD(1)=min(arr_minRDD);
    minRAD(1)=min(arr_minRAD);
end
%Find critical distances for the attackers
for i=1:NA-1
    for l=i+1:NA
        RAiAl(l)=norm(rA(:,i)-rA(:,l));
    end
    arr_minRAA(i)=min(RAiAl(i+1:end));
    %distance  of the attacker from the obstacles
    for kk=1:NO
        n=nO(kk);
        a=aO(kk);
        b=bO(kk);
        % RDjOl(l)=norm(rD(:,j)-rO(1:2,l));
        EAOl(kk)=(abs((rA(1,i)-rCO2(1,kk))/a))^(2*n)+(abs((rA(2,i)-rCO2(2,kk))/b))^(2*n)-1;
        EAO_crit(i,kk)=R_m_AO(i,kk)/EAOl(kk);
    end
end
minRAA(1)=min(arr_minRAA);
minEAO(1)=max(max(EAO_crit));
%%
dpsiAi0=0.02;
vA_ref_max=4;
global dpsi_var
syms dpsi_var
flagDefReachOpen=0;
flagDefReachClosed=0;
flagDefConnect=0;
flagAttInSight=0;
defReachCount=zeros(ND,1);
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

flagGather=1;
flagSeek=0;
flagEnclose=0;
flagHerd=0;

%%
%Initialize a ROS node to communicate with Gazebo
rosshutdown
rosinit('http://192.168.1.119:11311');

flagUseGroundTruth=1;

%Create subcribers and publishers for the attackers
for i=1:NA
    %Subcribe the odometry data for the attackers
    if flagUseGroundTruth==1
    formatspec1 = '/attacker%d/station/odometry';     
    formatspec2 = '/attacker%d/station/pose';
    else
        formatspec1 = '/attacker%d/odometry_sensor1/odometry';     
    formatspec2 = '/attacker%d/odometry_sensor1/pose';     
    end
    %formatspec = '/pelican/odometry_sensor1/odometry';
    string = sprintf(formatspec1,i);
    subA_odometry(i) = rossubscriber(string,'nav_msgs/Odometry');
   
    %formatspec = '/pelican/odometry_sensor1/pose';
    string = sprintf(formatspec2,i);
    subA_pose(i) = rossubscriber(string,'geometry_msgs/Pose');
    
    %Create a publisher for the raw control commands for the attackers
   formatspec3 = '/attacker%d/command/motor_speed';
    string = sprintf(formatspec3,i);
    pubA_rot_speed(i)= rospublisher(string,'mav_msgs/Actuators');
    msgPubA_rot_speed(i) = rosmessage(pubA_rot_speed(i));
    
    %Create a publisher for the pose of the attackers
    formatspec4 = '/attacker%d/command/pose';
    %formatspec = '/pelican/odometry_sensor1/pose';
    string = sprintf(formatspec4,i);
    pubA_pose (i)= rospublisher(string,'geometry_msgs/PoseStamped');
    msgPubA_pose(i)=rosmessage(pubA_pose(i));
end

%Create subcribers and publishers for the defenders
for j=1:ND
    %Subcribe the odometry data for the defenders
     if flagUseGroundTruth==1
    formatspec1 = '/defender%d/station/odometry';     
    formatspec2 = '/defender%d/station/pose';
    else
        formatspec1 = '/defender%d/odometry_sensor1/odometry';     
    formatspec2 = '/defender%d/odometry_sensor1/pose';     
    end
    
    %formatspec = '/pelican/odometry_sensor1/odometry';
    string = sprintf(formatspec1,j);
    subD_odometry(j) = rossubscriber(string,'nav_msgs/Odometry');
    string = sprintf( '/defender%d/station/odometry',j)
    subD_odometry_ground(j) = rossubscriber(string,'nav_msgs/Odometry');
     
    %formatspec = '/pelican/odometry_sensor1/pose';
    string = sprintf(formatspec2,j);
    subD_pose(j) = rossubscriber(string,'geometry_msgs/Pose');
    
    %Create a publisher for the raw control commands for the defebders
    formatspec3 = '/defender%d/command/motor_speed';
    string = sprintf(formatspec3,j);
    pubD_rot_speed(j)= rospublisher(string2,'mav_msgs/Actuators');
    msgPubD_rot_speed(j) = rosmessage(pubD_rot_speed(j));
    
    %Create a publisher for the pose of the defebders
    formatspec4 = '/defender%d/command/pose';
    %formatspec = '/pelican/odometry_sensor1/pose';
    string = sprintf(formatspec4,j);
    pubD_pose (j)= rospublisher(string,'geometry_msgs/PoseStamped');
    msgPubD_pose(j)=rosmessage(pubD_pose(j));
end

client_gazebo_pause_physics = rossvcclient('/gazebo/pause_physics');
client_gazebo_unpause_physics = rossvcclient('/gazebo/unpause_physics');
%rosTime=rossubscriber('/clock','rosgraph_msgs/Clock');


datum_time=rostime('now');
%call(client_gazebo_pause_physics);  %Pause Gazebo physics engine

flagSendMotorCommands=0;   %If 1 then run the launch file that does not lanuch the onboard controller

%Make the quadrotors hover at an altitude z=z_h
z_h=1.5;

dT=0.04; %Time step size for hovering controller
dt_i=0.01;  %time step size for the inner loop of the low level controller
T=0;
tcount=1;
tol=1;
flagReachedAltitude=0;
while (~flagReachedAltitude)
    %Subscribe to the data from the Gazebo
    for i=1:NA
        msgSubA_pose=subA_pose(i).LatestMessage;
        pos(:,i)=[msgSubA_pose.Position.X, msgSubA_pose.Position.Y, msgSubA_pose.Position.Z]';
        orientation(:,i)=[msgSubA_pose.Orientation.X, msgSubA_pose.Orientation.Y, msgSubA_pose.Orientation.Z, msgSubA_pose.Orientation.W]';
        
        msgSubA_odometry=subA_odometry(i).LatestMessage;
        Linear= msgSubA_odometry.Twist.Twist.Linear;
        vel(:,i)=[Linear.X,Linear.Y,Linear.Z]';
        Angular=msgSubA_odometry.Twist.Twist.Angular;
        ang_vel(:,i)=[Angular.X,Angular.Y,Angular.Z]';
    end
    
    for j=1:ND
        jj=j+NA;
        msgSubD_pose=subD_pose(j).LatestMessage;
        pos(:,jj)=[msgSubD_pose.Position.X, msgSubD_pose.Position.Y, msgSubD_pose.Position.Z]';
        orientation(:,jj)=[msgSubD_pose.Orientation.X, msgSubD_pose.Orientation.Y, msgSubD_pose.Orientation.Z, msgSubD_pose.Orientation.W]';
        
        msgSubD_odometry=subD_odometry(j).LatestMessage;
        Linear= msgSubA_odometry.Twist.Twist.Linear;
        vel(:,jj)=[Linear.X,Linear.Y,Linear.Z]';
        Angular=msgSubD_odometry.Twist.Twist.Angular;
        ang_vel(:,jj)=[Angular.X,Angular.Y,Angular.Z]';
    end
    
    %Desired quantities for the agents
    T=T+dT;
    acc_des=[zeros(3,N)];
    pos_des=[YA(1:2,:),XD(1:2,1:ND); 7+min(2*T,z_h-7)*ones(1,N)];
    vel_des=[YA(3:4,:),XD(3:4,1:ND); zeros(1,N)];
    yaw_des=[atan2(XA(2,2)-XA(2,1),XA(1,2)-XA(1,1)),zeros(1,N -1)];
    %yaw_des=zeros(1,N);
    yaw_rate_des=zeros(1,N);
    
    %Find the raw motor commands for the agents
    rotor_speeds=quad_controller(pos,vel,orientation,ang_vel,pos_des,vel_des,acc_des,yaw_des,yaw_rate_des);
    
    %Send commands to Gazebo    
    %For attackers
    if flagSendMotorCommands
        %Raw motor commands
        for i=1:NA
            msgPubA_rot_speed(i).AngularVelocities=rotor_speeds(:,i);
            send(pubA_rot_speed(i),  msgPubA_rot_speed(i)) ;
        end
    else
        %Position commands
        for i=1:NA
            msgPubA_pose(i).Pose.Position.X=XA(1,i);
            msgPubA_pose(i).Pose.Position.Y=XA(2,i);
            msgPubA_pose(i).Pose.Position.Z=z_h;
            if i==1
            rotMat=[cos(yaw_des(i)), sin(yaw_des(i)),0;  -sin(yaw_des(i)), cos(yaw_des(i)),0; 0,0,1]';
                quat=rotm2quat(rotMat);
                msgPubA_pose(i).Pose.Orientation.X=quat(2);
                msgPubA_pose(i).Pose.Orientation.Y=quat(3);
                msgPubA_pose(i).Pose.Orientation.Z=quat(4);
                msgPubA_pose(i).Pose.Orientation.W=quat(1);
            end
            send(pubA_pose(i),  msgPubA_pose(i)) ;
        end
    end   
    %For defenders    
    if flagSendMotorCommands
        %Raw motor commands
        for j=1:ND
            jj=j+NA;
            msgPubD_rot_speed(i).AngularVelocities=rotor_speeds(:,jj);
            send(pubD_rot_speed(j),  msgPubD_rot_speed(j)) ;
        end
    else
        for j=1:ND
            jj=j+NA;
            %Position commands
            msgPubD_pose(j).Pose.Position.X=XD(1,j);
            msgPubD_pose(j).Pose.Position.Y=XD(2,j);
            msgPubD_pose(j).Pose.Position.Z=z_h;
            send(pubD_pose(j),  msgPubD_pose(j)) ;
        end
    end
    %save the data
    hover_pos_history(:,tcount)=pos(:);
    hover_vel_history(:,tcount)=vel(:);
    hover_rotor_speeds_history(:,tcount)=rotor_speeds(:);
    tcount=tcount+1;
    
    %(Ensure that not much MATLAB computation is happening between unpausing and pausing to avoid unstable system behaviour)
    call(client_gazebo_unpause_physics);   %Unpause Gazebo physics engine
    pause(dT);
    %elapsedTime = rostime('now')-datum_time;
    %elapsedTime = rostime('now')-datum_time;
    %elapsedTime = rostime('now')-datum_time;
    call(client_gazebo_pause_physics);  %Pause Gazebo physics engine
    %    dT=double(elapsedTime.Sec)+double(elapsedTime.Nsec)*10^(-9);
    
    agentCount=0;
    for i=1:N
        if norm(pos(:,i)-pos_des(:,i))<tol %&& norm(orientation(:,i)-[0,0,0,1]')<0.2
            agentCount=agentCount+1;
        end
    end
    if agentCount>N-1
        flagReachedAltitude=1;
    end
end

yaw_des0=yaw_des;

for ti=1:Niter
    if mod(ti,1000)==0
        ti
    end
    
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
    countDInS=0;
    for jj=1:ND
        if norm(rD(:,jj)-rS)<rho_S
            countDInS=countDInS+1;
        end
    end
    if countAInS>NA-1 && countDInS>ND-1
        break;
    end
    if safe_flag==1
        if t-t_safe<=Delta_Tt
            psi=psi+(t-t_safe)*(pi/2-epsilon)/(Delta_Tt);
            psi0=psi;
            t0=t;
            %rArS0=norm(rA(:)-rS);
            %psi_dot=psi_dot+(pi/2-epsilon)/(Delta_Tt);
        else
            psi=psi+pi/2-epsilon;
        end
    end
    
    %Initialize the desired position to the current position
    xd0=XD(:,1:ND);
    XD_Des(:,ti+1)=xd0(:);
    XD_Des_dot(:,ti+1)=zeros(ND*4,1);
    
    pairDA=[1:ND]';   %pairDA=ones(ND,1);
    
    kref=0.1;
    
    refTraj=ones(4,ND);
    RefTraj(:,ti)=refTraj(:);
    
    %%
    %desired positions for the leader attackers
    XA_lead_goal=[rP;0;0];
    XA_goal(:,1)=XA_lead_goal;
    for i=2:NA
        thetaA=360-360*(i-1)/NA;
        
        XA_goal(1:2,i)=XA(1:2,1)+rA_follow(1:2,i);
        XA_goal(3:4,i)=XA(3:4,1);
        XA_goal_dot(:,i)=[XA_goal(3:4,i); uA(1:2,1)];
        if flagEnclose==1 || flagHerd==1
            if i>4
                XA_goal(:,i)=XA_lead_goal;
                XA_goal_dot(:,i)=zeros(4,1);
            end
        end
        %For flocking
        % XA_goal(1:2,i)=rP;
    end
    
    %[uA, uA0, R_AO_min, R_AAProjS_min, vA_des, vA_des_dot, rAProj1, flagAttInObs, betaAv0, speedA0, mAv0, cAv0, SigmaProdD, FA, FA_dot] = controlAttacker3(XA,XA_goal,vA,flagAttInObs,betaAv0,speedA0,mAv0,cAv0,XD,WA,WDString,NA,ND);
    
    %When the attackers are to be split
    [uA, uA0, R_AO_min, R_AAProjS_min, vA_des, vA_des_dot, rAProj1, flagAttInObs, betaAv0, speedA0, mAv0, cAv0, SigmaProdD, FA, FA_dot] = controlAttacker4(XA,XA_goal,XA_goal_dot,vA,flagAttInObs,betaAv0,speedA0,mAv0,cAv0,XD,WA,WDString,NA,ND);
    rAProj_arr(:,:,ti)=rAProj1;  %row=attacker, column=obstacle, slice=time
    SigmaProdD_arr(:,ti)=SigmaProdD;
    
    
    %Control action for the defenders
    %[uD, arr_EDO,e_vD_des]= controlDefender4(XD,indDef, XD_des, XD_des_dot, XA,pairDA,refTraj, ND);
    rAcm=sum(XA(1:2,:),2)/NA;
    vAcm=sum(XA(3:4,:),2)/NA;
    rAcm_arr(:,ti)=rAcm;
    vAcm_arr(:,ti)=vAcm;
    
    %Find if the defenders are connected with each other
    WDString=zeros(ND);
    for j=1:ND-1
        j1=find(assign==j);
        j2=find(assign==j+1);
        if norm(XD(1:2,j1)-XD(1:2,j2))<=R_DD_string
            WDString(j1,j2)=1;
            WDString(j2,j1)=1;
        else
            WDString(j1,j2)=0;
            WDString(j2,j1)=0;
        end
    end
    if(0)
        j1=find(assign==1); %Defender corresponding to first desired position
        jND=find(assign==ND);   %Defender corresponding to last desired position
        if flagGather~=1 && flagSeek~=1
            if norm(XD(1:2,j1)-XD(1:2,jND))<=R_DD_string
                WDString(j1,jND)=1;
                WDString(jND,j1)=1;
            else
                WDString(j1,jND)=0;
                WDString(jND,j1)=0;
            end
        end
    end
    
    if flagGather==1 && flagSeek~=1 && flagEnclose~=1 && flagHerd~=1
        % ~(norm(rAcm-rDcm)<(2*rho_sn) && flagAttInSight) && flagHerd~=1 && flagDForm~=1 %&& norm(vAcm)>0.1   %%((rAcm-rDcm)'*[cos(thetaAcm0);sin(thetaAcm0)])<0
        XD_des=XD_des0;
        XD_des_dot=XD_des_dot0;
        
        %Time optimal control
        [uD, arr_EDO] = controlDefender5(XD, SD, [1:ND+1], assign, XD_des, XD_des_dot, assignedPath, assignedPathVel, time(ti), startTime, ND );
        ti_1=ti;
        
        for j=1:ND
            if norm(XD(1:2,indDef(j))-XD_des(1:2,j))<1 %&& norm(XD(3:4,j))<1e-8 && flagDefReachOpen~=1
                defReachCount(j)=1;
                if sum(defReachCount)>=ND
                    flagDefReachOpen=1;
                    flagSeek=1;
                    flagGather=0;
                    break;
                end
            end
        end
        ti_g=ti;
    elseif flagGather~=1 && flagSeek==1 && flagEnclose~=1 && flagHerd~=1
        
        j1=find(assign==1); %Defender corresponding to first desired position
        jND=find(assign==ND);
        %[uD, arr_EDO,FlagDefInObs,betaDv0,mDv0,cDv0]= controlDefender3(XD,indDef, XD_des, XD_des_dot,FlagDefInObs,betaDv0,mDv0,cDv0, XA,pairDA,refTraj, ND);
        
        [XD_des,XD_des_dot,phi_dot,phi_ddot,uDFc_trans,flagAttInSight]=defDesiredOpenForm(XD(:,ND+1),RDF_open,YA,phi,phi_dot,NA,ND);
        phi=phi+dt*phi_dot;
        if phi>2*pi
            phi=phi-2*pi;
        end
        phi_dot=phi_dot+dt*phi_ddot;
        [uD, arr_EDO] = controlFiniteTimeTrajTracking(XD,indDef, XD_des, XD_des_dot, uDFc_trans, YA, ND, 1 );

        %uD(:,ND+1)=uDFc_trans;
        %uD=zeros(2,ND+1);
        ti_2=ti;
        
        if (norm(rAcm-rDcm)<(2*rho_sn) && flagAttInSight) %&& flagHerd~=1 && flagDForm~=1 %&& norm(vAcm)>0.1
            flagEnclose=1;
            flagSeek=0;
        end
        ti_s=ti;
    elseif flagGather~=1 && flagSeek~=1 && flagEnclose==1 && flagHerd~=1 %flagHerd~=1  && flagDefConnect~=1 %countDDes<ND  %Enclosing phase
        %thetaAcm=atan2(vAcm(2),vAcm(1));
        if flagDForm~=1
            flagDForm=1;
        end
        for j=1:ND
            RD0=RDF_closed;%3*rDmin;
            %phi=desPhiForClosedForm(rAcm);
            if flagNDeven~=1
                thetaD=phi+2*pi*(j)/ND-pi/ND;
            else
               thetaD=phi+2*pi*(j)/ND-pi/ND;
            end
            rD_des(:,j)=rAcm+RD0*[cos(thetaD), sin(thetaD);]';
            XD_des(:,j)=[rD_des(:,j);vAcm];
            XD_des_dot(:,j)=[0,0,0,0]';
            rSD_goal(:,j)=rS+RD0*[cos(thetaD), sin(thetaD);]';
        end
        
        %Check if the terminal defenders should be connected or not
        j1=find(assign==1); %Defender corresponding to first desired position
        jND=find(assign==ND);   %Defender corresponding to last desired position
        if (norm(XD(1:2,j1)-XD_des(1:2,1))< bd && norm(XD(1:2,jND)-XD_des(1:2,ND))<bd) ||  (norm(XD(1:2,j1)-XD(1:2,jND))<=0.97*R_DD_string && inhull(rAcm',XD(1:2,1:ND)'))
            WDString(j1,jND)=1;
            WDString(jND,j1)=1;
        else
            WDString(j1,jND)=0;
            WDString(jND,j1)=0;
        end
        
        %Check if all the defenders are connected to each other and the
        %StringNet is formed
        countDefConnect=0;
        for j=1:ND
            if j>=ND
                if WDString(indDef(ND),indDef(1))==1
                    countDefConnect=countDefConnect+1;
                end
            else
                if WDString(indDef(j),indDef(j+1))==1
                    countDefConnect=countDefConnect+1;
                end
            end
        end
        if countDefConnect==ND
            flagDefConnect=1;
            flagHerd=1;
            flagEnclose=0;
            rDcm=sum(XD(1:2,:),2)/ND;
            thetaD1=atan2(XD(2,indDef(1))-rDcm(2),XD(1,indDef(1))-rDcm(1));
            for j=1:ND
                thetaD=thetaD1+2*pi*(j-1)/(ND);
                XD_herd_des0(:,indDef(j))=rDcm+Rjj_bar_closed(1,ND+1)*[cos(thetaD),sin(thetaD)]';
            end
        end
        %%Use the assignment from the previous phase to continue
        %%completing the StringNet
        %                  XD_des=XD_des(:,assign);   %change the XD_des according to the assignment
        %                  XD_des_dot=XD_des_dot(:,assign);
        
        [ uD, arr_EDO] = controlDefenderFormation4(XD, indDef, assign, XD_des, XD_des_dot, uDFc_trans, YA,NA,ND,1);
        % uD=[uD,[0;0]]; %For the virtual defender at the CoM of the attackers
        % indDef(assign)=indDef;   %Indices of the defenders corresponding to the formation in that order
        
        
        %  [uD, arr_EDO] = controlDefenderFormation2(XD, indDef, XD, WD_closed, WDString, Rjj_bar_closed, XA,NA,ND);
        ti_3=ti;
        ti_e=ti;
        %XDFc=[rAcm;vAcm];
        XDF_des=XD_des;
        X(4*(NA+ND+1)-3:4*(NA+ND+1),ti)=[rAcm;vAcm]; %The virtual defender's state
    elseif flagDefReachClosed~=1
        
        %Check if the terminal defenders should be connected or not
        j1=find(assign==1); %Defender corresponding to first desired position
        jND=find(assign==ND);   %Defender corresponding to last desired position
        if norm(XD(1:2,j1)-XD(1:2,jND))<=R_DD_string
            WDString(j1,jND)=1;
            WDString(jND,j1)=1;
        else
            WDString(j1,jND)=0;
            WDString(jND,j1)=0;
        end
        
        flagDefReachClosed=1;
        flagHerd=1;
        for j=1:ND
            %             if norm(XD(1:2,indDef(j))-XD_herd_des0(1:2,j))<1e-3 && norm(XD(3:4,j))<1e-8 && flagDefReachClosed~=1
            if norm(XD(1:2,indDef(j))-XDF_des(1:2,j))<2 && norm(XD(3:4,j))<1e-1 && flagDefReachClosed~=1
                defReachCount(j)=1;
                if sum(defReachCount)>=ND
                    flagDefReachClosed=1;
                    flagHerd=1;
                    break;
                end
            end
        end
        
        uD=zeros(2,ND+1);
        for j=1:ND
            uD(:,indDef(j))=-0.1*(XD(1:2,indDef(j))-XDF_des(1:2,j))-0.1*(XD(3:4,indDef(j)));
            infNorm_uD=max(abs(uD(1,indDef(j))),abs(uD(2,indDef(j))));
            if infNorm_uD>u_maxD(indDef(j))
                uD(:,indDef(j))=uD(:,indDef(j))*u_maxD(indDef(j))/infNorm_uD;
            end
        end
        uD(:,ND+1)=C_d*XD(3:4,ND+1);
        %uD(:,ND+1)=sum(uD(:,1:ND),2)/ND;
        %XD(:,ND+1)=sum(XDF_des,2)/ND;
        XD_des=XDF_des;
        ti_3=ti;
        ti_e=ti;
    elseif flagGather~=1 && flagSeek~=1 && flagEnclose~=1 && flagHerd==1
        %Check if the terminal defenders should be connected or not
        j1=find(assign==1); %Defender corresponding to first desired position
        jND=find(assign==ND);   %Defender corresponding to last desired position
        if norm(XD(1:2,j1)-XD(1:2,jND))<=R_DD_string
            WDString(j1,jND)=1;
            WDString(jND,j1)=1;
        else
            WDString(j1,jND)=0;
            WDString(jND,j1)=0;
        end
        
        %[uD, arr_EDO]= controlDefender4(XD,indDef, XD_des, XD_des_dot, XA, WDString, Rjj_bar_closed,ND);
        
        %Find the desired positions of the defenders on the closed
        %formation that moves towards the safe area
        if ti-ti_3<100
            XD(:,ND+1)=[rAcm;vAcm];
        end
        XDFc_des=XD(:,ND+1);
        delta_t=dt*(ti-ti_3);
        [XD_des,XD_des_dot,uDFc_trans]=defDesiredClosedForm(XDFc_des,RDF_closed,phi,YA,NA,ND,flagNDeven,delta_t);
        
        [ uD, arr_EDO] = controlFiniteTimeTrajTracking(XD,indDef, XD_des, XD_des_dot, uDFc_trans, YA, ND, 0 );
        %uD(:,ND+1)=XD_des_dot(3:4,1);
        
    end
    
    XD_Des(:,ti+1)=XD_des(:);
    WDString_mat(:,:,ti)=WDString;
    
     U(:,ti)=[uA(:);uD(:)];
     X(:,ti+1)=modifiedDIDynamics(X(:,ti),U(:,ti));
    
    %Run the inner loop at more frequency than the outer one.
    dt_i0=0;
    while dt_i0<dt
        %Desired quantities for the agents
        %acc_des=[uA,uD(:,1:ND); zeros(1,N)];
        Xti1=reshape(X(:,ti+1),4,N+1);
        acc_des=[uA-C_d*Xti1(3:4,1:NA),uD(:,1:ND)-C_d*Xti1(3:4,NA+1:N); zeros(1,N)];
        pos_des=[Xti1(1:2,1:N); z_h*ones(1,N)];
        vel_des=[Xti1(3:4,1:N); zeros(1,N)];
        yaw_des=[atan2(XA(2,2)-XA(2,1),XA(1,2)-XA(1,1)),zeros(1,NA -1),0, atan2(vel(2,NA+2),vel(1,NA+2)), zeros(1,ND-2)];
        %yaw_des=yaw_des0;
        yaw_rate_des=zeros(1,N);
        acc_des_history(:,ti)=acc_des(:);
        %Find the raw motor commands for the agents
        %rotor_speeds=quad_controller(pos,vel,orientation,ang_vel,pos_des,vel_des,acc_des,yaw_des,yaw_rate_des);
        
        %Send commands to Gazebo
        %Commands for attackers        
        if flagSendMotorCommands
            %Raw motor commands
            for i=1:NA
                msgPubA_rot_speed(i).AngularVelocities=rotor_speeds(:,i);
                send(pubA_rot_speed(i),  msgPubA_rot_speed(i)) ;
            end
        else
            %Position commands
            for i=1:NA                
                msgPubA_pose(i).Pose.Position.X=X(4*i-3,ti+1);
                msgPubA_pose(i).Pose.Position.Y=X(4*i-2,ti+1);
                msgPubA_pose(i).Pose.Position.Z=z_h;
                rotMat=[cos(yaw_des(i)), sin(yaw_des(i)),0;  -sin(yaw_des(i)), cos(yaw_des(i)),0; 0,0,1]';
                quat=rotm2quat(rotMat);
                msgPubA_pose(i).Pose.Orientation.X=quat(2);
                msgPubA_pose(i).Pose.Orientation.Y=quat(3);
                msgPubA_pose(i).Pose.Orientation.Z=quat(4);
                msgPubA_pose(i).Pose.Orientation.W=quat(1);
                send(pubA_pose(i),  msgPubA_pose(i)) ;
            end
        end
        %Commands for defenders       
        if flagSendMotorCommands
             %Raw motor commands
            for j=1:ND
                jj=j+NA;
                msgPubD_rot_speed(i).AngularVelocities=rotor_speeds(:,jj);
                send(pubD_rot_speed(j),  msgPubD_rot_speed(j)) ;
            end
        else
            %Position commands
            for j=1:ND
                jj=j+NA;                
                msgPubD_pose(j).Pose.Position.X=X(4*jj-3,ti+1);
                msgPubD_pose(j).Pose.Position.Y=X(4*jj-2,ti+1);
                msgPubD_pose(j).Pose.Position.Z=z_h;
                send(pubD_pose(j),  msgPubD_pose(j)) ;
            end
        end
        call(client_gazebo_unpause_physics);   %Unpause Gazebo physics engine
        %elapsedTime = rostime('now')-datum_time;
        pause(0.04);
        call(client_gazebo_pause_physics);  %Pause Gazebo physics engine
        
        %Subscribe to the data from the Gazebo
        for i=1:NA
            msgSubA_pose=subA_pose(i).LatestMessage;
            pos(:,i)=[msgSubA_pose.Position.X, msgSubA_pose.Position.Y, msgSubA_pose.Position.Z]';
            orientation(:,i)=[msgSubA_pose.Orientation.X, msgSubA_pose.Orientation.Y, msgSubA_pose.Orientation.Z, msgSubA_pose.Orientation.W]';
            
            msgSubA_odometry=subA_odometry(i).LatestMessage;
            Linear= msgSubA_odometry.Twist.Twist.Linear;
            vel(:,i)=[Linear.X,Linear.Y,Linear.Z]';
            Angular=msgSubA_odometry.Twist.Twist.Angular;
            ang_vel(:,i)=[Angular.X,Angular.Y,Angular.Z]';
        end
        %Measurements of the attackers
        YA=[pos(1:2,1:NA);vel(1:2,1:NA)];
        rA=YA(1:2,:);
        vA=YA(3:4,:);
        
        for j=1:ND
            jj=j+NA;
            msgSubD_pose=subD_pose(j).LatestMessage;
            pos(:,jj)=[msgSubD_pose.Position.X, msgSubD_pose.Position.Y, msgSubD_pose.Position.Z]';
            orientation(:,jj)=[msgSubD_pose.Orientation.X, msgSubD_pose.Orientation.Y, msgSubD_pose.Orientation.Z, msgSubD_pose.Orientation.W]';
            
            msgSubD_odometry=subD_odometry(j).LatestMessage;
            Linear= msgSubA_odometry.Twist.Twist.Linear;
            vel(:,jj)=[Linear.X,Linear.Y,Linear.Z]';
            Angular=msgSubD_odometry.Twist.Twist.Angular;
            ang_vel(:,jj)=[Angular.X,Angular.Y,Angular.Z]';
        end
        XDp=XD;    %previous state
        XD=[[pos(1:2,NA+1:end);vel(1:2,NA+1:end)],X(4*N+1:4*N+4,ti+1)];
        
        rD=XD(1:2,:);
        vD=XD(3:4,:);
        dt_i0=dt_i0+dt_i;
    end
    
    pos_history(:,ti+1)=pos(:);
    vel_history(:,ti+1)=vel(:);
    rotor_speeds_history(:,ti+1)=rotor_speeds(:);
    
    %Saturate the velocity if beyond the maximum
    for i=1:NA
        norm_vA=norm(vA(:,i));
        if norm_vA>v_maxA(i)
            vA(:,i)=vA(:,i)*v_maxA(i)/norm_vA;
        end
    end
    for j=1:ND
        norm_vD=norm(vD(:,j));
        if norm_vD>v_maxD(j)
            vD(:,j)=vD(:,j)*v_maxD(j)/norm_vD;
        end
    end
    XA=[rA;vA];
    XD=[rD;vD];
    %X(:,ti+1)=[XA(:);XD(:)];
    
    XGazebo(:,ti+1)=[XA(:);XD(:)];
    
    rDcm=sum(XD(1:2,1:ND),2)/ND;
    vDcm=sum(XD(3:4,1:ND),2)/ND;
    
    %Distance travelled (position along the path)
    for j=1:ND
        SD(j)=SD_arr(j,ti)+norm(XD(1:2,j)-XDp(1:2,j));
    end
    SD_arr(:,ti+1)=SD;
    %X(4*NA+1:4*N,ti+1)=XD(:);
    %%
    %         XD=XD_des;
    %          X(4*NA+1:end,ti+1)=XD(:);
    %if ND>1
    arr_minRDD=[];
    arr_minRAD=[];
    arr_minRAA=[];
    for j=1:ND
        if j<ND
            for l=j+1:ND
                RDjDl(l)=norm(rD(:,j)-rD(:,l));
            end
            arr_minRDD(j)=min(RDjDl(j+1:end));
        end
        for i=1:NA
            arr_minRAD(i,j)=norm(rD(:,j)-rA(:,i));
        end
        %relative distance of the defenders from Ok
        for k=1:NO
            if(0)
                E_DOk=Inf;
                for kk=1:length(rVO{k}(1,:))-1
                    dx=rVO{k}(1,kk+1)-rVO{k}(1,kk);
                    dy=rVO{k}(2,kk+1)-rVO{k}(2,kk);
                    mDD=dy/dx;  %slope of line joining the defenders
                    cDD=rVO{k}(2,kk)-mDD*rVO{k}(1,kk);
                    rDProjO(1,1)=(mDD*rD(2,j)+rD(1,j)-mDD*cDD)/(1+mDD^2);
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
                    E_DOk=min(E_DOk,norm(rD(:,j)-rDProjO));
                end
                E_DOk=E_DOk-rho_safe;
                arr_EDO(j,k)=E_DOk;
            end
            arr_EDO_crit(j,k)=R_m_DD/arr_EDO(j,k);
        end
    end
    if ND>1
        minEDO(ti)=max(max(arr_EDO_crit));
    end
    if ~isempty(arr_minRDD)
        minRDD(ti+1)=min(arr_minRDD);
        
    end
    
    for i=1:NA
        arr_minRAD(i,ND)=norm(rD(:,ND)-rA(:,i));
    end
    minRAD(ti+1)=min(min(arr_minRAD));
    
    %Find critical distances for the attackers
    for i=1:NA
        if i<NA
            for ii=i+1:NA
                RAiAl(ii)=norm(rA(:,i)-rA(:,ii));
            end
            arr_minRAA(i)=min(RAiAl(i+1:end));
        end
    end
    if ~isempty(arr_minRAA)
        minRAA(ti+1)=min(arr_minRAA);
    end
    if NA>1
        minEAO(ti)=rhoAD_safe/(R_AO_min+rhoAD_safe);
    end
    minRAAProjS(ti)=R_m_AD/R_AAProjS_min;
    
    t=t+dt;
    time(ti+1)=t;
    %ti=ti+1;
end
if ti==Niter
    ti=ti+1
end
%%
vA_dot_arr(:,ti)=vA_dot_arr(:,ti-1);
FA_dot_arr(:,ti)=FA_dot_arr(:,ti-1);
U(:,ti)=U(:,ti-1);
%vA_Des(:,ti)=vA_Des(:,ti-1);
FA_arr(:,ti)=FA_arr(:,ti-1);
RefTraj(:,ti)=RefTraj(:,ti-1);
rAProj_arr(:,:,ti)=rAProj_arr(:,:,ti-1);
SigmaProdD_arr(:,ti)=SigmaProdD;
%e_vD_des_Arr(:,ti)=e_vD_des;
% Psi_prime(ti)=Psi_prime(ti-1);
% Psi_prime_dot(ti)=Psi_prime_dot(ti-1);
% Psi_prime_ddot(ti)=Psi_prime_ddot(ti-1);
% Psi(ti)=Psi(ti-1);
% Psi_des0(ti)=Psi_des0(ti-1);
% Psi_des(ti)=Psi_des(ti-1);
% Psi_des_dot(ti)=Psi_des_dot(ti-1);
minEDO(ti)=minEDO(ti-1); %not so true but won't make much difference for now, maybe explicitly calculate it, if required.
minEAO(ti)= minEAO(ti-1);
minRAAProjS(ti)=minRAAProjS(ti-1);
if (1) %delete the unnecessry elements ti<Niter
    time(:,ti+1:end)=[];
    X(:,ti+1:end)=[];
    U(:,ti+1:end)=[];
    XD_Des(:,ti+1:end)=[];
    vA_dot_arr(:,ti+1:end)=[];
    FA_dot_arr(:,ti+1:end)=[];
    minRDD(:,ti+1:end)=[];
    minRAD(:,ti+1:end)=[];
    minRAA(:,ti+1:end)=[];
    minEDO(:,ti+1:end)=[];
    minEAO(:,ti+1:end)=[];
    SigmaProdD_arr(:,ti+1:end)=[];
    rAProj_arr(:,:,ti+1:end)=[];
end
%minEAO(ti)=minEAO(ti-1);
%%

%Plot the trajectories
fontSize=18;
plotTraj(X(:,1:ti),rP,rS,rVO,rho_S,time(1:ti),NA,ND,'o','--',1,1,1,1);
set(gca,'fontsize',fontSize)
%%
%animateTraj(X,rP,rS,rVO,rho_Acon,rho_S,time,NA,ND,XD_Des,WDString_mat);
%%
for i=1:NA+ND
    colors{i}=rand(1,3);
end

plotControls=1;
if plotControls==1
    %Controls for attackers
    if(0)
        figure
        for j=1:N
            if j<=NA
                plot(time,U(2*(j-1)+1,:),'r-');
                hold on;
            else
                plot(time,U(2*(j-1)+1,:),'b-');
                hold on;
            end
            hold on;
        end
        xlabel('t [s]');
        ylabel('$u_x$')
        %Controls for defenders
        figure
        for j=1:N
            if j<=NA
                plot(time,U(2*(j-1)+2,:),'r-');
                hold on;
            else
                plot(time,U(2*(j-1)+2,:),'b-');
                hold on;
            end
            hold on;
        end
        xlabel('t [s]');
        ylabel('$u_y$')
    end
end

%norm of the control for all agents
if(1)
    figure('units','normalized','outerposition',[.2 0.4 .5 .4])
    
    sp2=subplot(2,1,2)
    sp2.Position(4)=0.07;
    hold on
    fill(dt*[0,ti_g,ti_g,0,0],u_maxD(1)*[1.1,1.1,1.3,1.3,1.1],'c','facealpha',0.9)
    fill(dt*[ti_g+1,ti_s,ti_s,ti_g+1,ti_g+1],u_maxD(1)*[1.1,1.1,1.3,1.3,1.1],'k','facealpha',0.9)
    fill(dt*[ti_s+1,ti_e,ti_e,ti_s+1,ti_s+1],u_maxD(1)*[1.1,1.1,1.3,1.3,1.1],'m','facealpha',0.9)
    fill(dt*[ti_e+1,ti,ti,ti_e+1,ti_e+1],u_maxD(1)*[1.1,1.1,1.3,1.3,1.1],'g','facealpha',0.9)
    axis off
    
    sp1=subplot(2,1,1)
    sp1.Position(2)=0.07+0.3;
    sp1.Position(4)=0.5;
    for j=1:N
        plot(time(1:ti),sqrt(U(2*(j-1)+1,1:ti).^2+U(2*(j-1)+2,1:ti).^2),'color',colors{j});
        hold on;
        if j<NA+1
            legends{j}=['$\mathcal{A}_{',num2str(j),'}$'];
        else
            legends{j}=['$\mathcal{D}_{',num2str(j-NA),'}$'];
        end
    end
    plot(time(1,1:ti),u_maxD(1)*ones(1,ti),'color',[1,0,0],'linestyle','--')
    plot(time(1,1:ti),u_maxA(1)*ones(1,ti),'color',[0,0,1],'linestyle','--')
    %     text(dt*ti_g/10,1.4*u_maxD(1),'$Gathering$','fontsize',fontSize,'color','c')
    %     text(dt*ti_g/10,1.4*u_maxD(1),'$Gathering$','fontsize',fontSize,'color','c')
    
    xlabel('t [s]');
    ylabel('$||u||, [m/s]$')
    legends{N+1}=['$\bar{u}_{d}$'];
    legends{N+2}=['$\bar{u}_{a}$'];
    leg = legend(legends,'orientation','horizontal')
    leg.ItemTokenSize =[14,10];
    ylim([-1,u_maxD(1)+5])
    %axis([0,dt*ti+10,-0.2,1.1*u_maxD(1)])
    set(gca,'fontsize',fontSize)
    
    
    
    
    %Plot minimum distances
    figure('units','normalized','outerposition',[.2 0.4 .5 .4])
    
    sp2=subplot(2,1,2)
    sp2.Position(4)=0.07;
    hold on
    fill(dt*[0,ti_g,ti_g,0,0],[1.1,1.1,1.3,1.3,1.1],'c','facealpha',0.9)
    fill(dt*[ti_g+1,ti_s,ti_s,ti_g+1,ti_g+1],[1.1,1.1,1.3,1.3,1.1],'k','facealpha',0.9)
    fill(dt*[ti_s+1,ti_e,ti_e,ti_s+1,ti_s+1],[1.1,1.1,1.3,1.3,1.1],'m','facealpha',0.9)
    fill(dt*[ti_e+1,ti,ti,ti_e+1,ti_e+1],[1.1,1.1,1.3,1.3,1.1],'g','facealpha',0.9)
    axis off
    sp1=subplot(2,1,1)
    sp1.Position(2)=0.07+0.3;
    sp1.Position(4)=0.5;
    plot(time,R_m_DD./minRDD,'b')
    hold on
    plot(time,R_m_AD./minRAD,'k')
    hold on
    plot(time,R_m_AA./minRAA,'r')
    plot(time,minEDO,'c')
    hold on
    plot(time,minEAO,'m',time,ones(size(time)),'r--')
    
    ylabel('Distance Ratios')
    xlabel('t [s]')
    ylh = get(gca,'ylabel');
    gyl = get(ylh);                                                         % Object Information
    ylp = get(ylh, 'Position');
    set(ylh, 'Rotation',90, 'Position',[-18,0.1,1], 'VerticalAlignment','middle', 'HorizontalAlignment','left')
    %legend('$\displaystyle \max_{k \in I_o} \displaystyle \max_{j \in I_d} \frac{E_{ok}^{djm}}{E_{ok}^{dj}(\textbf{r}_{dj})}$', '$\displaystyle \max_{k \in I_o}\frac{E_{ok}^{m}}{ E_{ok}(\textbf{r}_a)}$','y=1')
    leg=legend('$\mathcal{R}_{cr}^{dd}$','$\mathcal{R}_{cr}^{ad}$','$\mathcal{R}_{cr}^{aa}$','$\mathcal{R}_{cr}^{do}$','$\mathcal{R}_{cr}^{ao}$','y=1','orientation','horizontal')
    set(gca,'fontsize',fontSize)
    leg.ItemTokenSize =[14,10];
    ylim([0,1.2])
end

%%
if(0)
    %normal distances
    figure
    %subplot(2,1,1)
    plot(time,R_m_DD./minRDD,'b')
    %ylabel('$\displaystyle \max_{j\neq l \in I_d} \frac{R_{dd}^m}{||\textbf{r}_{dj}-\textbf{r}_{dl}||}$')
    %subplot(2,1,2)
    hold on
    plot(time,R_m_AD./minRAD,'r',time,ones(size(time)),'r--')
    %ylabel('$\displaystyle \max_{j \in I_d} \frac{R_{ad}^m}{||\textbf{r}_{dj}-\textbf{r}_{a}||}$')
    ylabel('Relative Distance')
    xlabel('t [s]')
    legend('$\displaystyle \max_{j\neq l \in I_d} \frac{R_{dd}^m}{||\textbf{r}_{dj}-\textbf{r}_{dl}||}$','$\displaystyle \max_{j \in I_d} \frac{R_{ad}^m}{||\textbf{r}_{dj}-\textbf{r}_{a}||}$','fontsize',44)
    
    %Elliptic distances
    figure
    %subplot(2,1,1)
    plot(time,minEDO,'b')
    %ylabel('$\displaystyle \min_{k \in I_o} \displaystyle \min_{j \in I_d} E_{ok}^{dj}(\textbf{r}_{dj})$')
    %subplot(2,1,2)
    hold on
    plot(time,minEAO,'r',time,ones(size(time)),'r--')
    %ylabel('$\displaystyle \min_{k \in I_o} E_{ok}(\textbf{r}_a)$')
    ylabel('Relative Elliptic Distance')
    xlabel('t [s]')
    legend('$\displaystyle \max_{k \in I_o} \displaystyle \max_{j \in I_d} \frac{E_{ok}^{djm}}{E_{ok}^{dj}(\textbf{r}_{dj})}$', '$\displaystyle \max_{k \in I_o}\frac{E_{ok}^{djm}}{ E_{ok}(\textbf{r}_a)}$','y=1')
end



%Plot defenders trajectory vs their desired trajectories
plotDesiredTraj=0;
if plotDesiredTraj==1
    nA=4*NA;
    for j=1:ND
        figure
        subplot(2,1,1)
        plot(time,XD_Des(4*j-3,:),time,X(nA+4*j-3,:),time,X(4*j-3,:))
        ylabel('$x_{d,err}$ [m]')
        %legend('$x_{des}$','$x$')
        subplot(2,1,2)
        plot(time,XD_Des(4*j-2,:),time,X(nA+4*j-2,:))
        ylabel('$y_{d,err}$ [m]')
        xlabel('t [s]')
        
        figure
        subplot(2,1,1)
        plot(time,XD_Des(4*j-1,:),time,X(nA+4*j-1,:))
        ylabel('$v_x [m/s]$')
        legend('$v_{x_{des}}$','$v_{x}$')
        subplot(2,1,2)
        plot(time,XD_Des(4*j,:),time,X(nA+4*j,:))
        ylabel('$v_y [m/s]$')
        xlabel('t [s]')
    end
end
plotSpeed=0;
if plotSpeed==1
    nA=4*NA;
    for j=1:ND
        figure
        % subplot(2,1,1)
        plot(time,sqrt(X(nA+4*j,:).^2+X(nA+4*j-1,:).^2))
        ylabel('$||v||$')
        xlabel('t [s]')
    end
end

%Plot attackers trajectory vs their desired trajectories
plotDesiredTrajA=0;
if plotDesiredTrajA==1
    for i=1:NA
        figure
        subplot(2,1,1)
        plot(time,X(4*i-3,:))
        ylabel('$x_{a}$ [m]')
        %legend('$x_{des}$','$x$')
        subplot(2,1,2)
        plot(time,X(4*i-2,:))
        ylabel('$y_{a}$ [m]')
        xlabel('t [s]')
        
        figure
        subplot(2,1,1)
        plot(time,vA_Des(2*i-1,:)-X(4*i-1,:))
        ylabel('$v_{xa_{err}}$ [m/s]')
        %legend('$x_{des}$','$x$')
        subplot(2,1,2)
        plot(time,vA_Des(2*i,:)-X(4*i,:))
        ylabel('$v_{ya_{err}}$ [m/s]')
        xlabel('t [s]')
    end
end



%Plot sigmaProdD
if(1)
    figure
    for i=1:NA
        plot(time,SigmaProdD_arr(i,:))
        ylabel('$sigmaProdD$')
        xlabel('t [s]')
        hold on;
    end
end
