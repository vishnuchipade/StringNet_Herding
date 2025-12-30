function [rotor_speeds]=quad_controller(pos,vel,orientation,ang_vel,pos_des,vel_des,acc_des,yaw_des,yaw_rate_des)
% ------------------------------
% Author: Vishnu S. Chipade
% Email: vishnuc@umich.edu
% ------------------------------


global params
m=params.m;
g=params.g;
I=params.I;
alloc_mat=params.alloc_mat;
inv_alloc_mat=params.inv_alloc_mat;
acc2omega=params.acc2omega;
pos_gain=params.pos_gain;
vel_gain=params.vel_gain;
att_gain=params.att_gain;
ang_rate_gain=params.ang_rate_gain;
att_gain_mar=params.att_gain_mar;
ang_rate_gain_mar=params.ang_rate_gain_mar;
omega_max=params.omega_max;

for i=1:length(pos(1,:))

quat=[orientation(4,i);orientation(1:3,i)];  %quaternion vector
        
rot_mat=quat2rotm(quat');
vel_W=rot_mat*vel;
pos_err=pos(:,i)-pos_des(:,i);
vel_err=vel_W(:,i)-vel_des(:,i);
acc_W=-[0,0,g]'+pos_gain*pos_err/m+vel_gain*vel_err/m-acc_des(:,i);

%Angular Moments
b1_des=[cos(yaw_des(i)), sin(yaw_des(i)),0]';
b3_des=-acc_W/(norm(acc_W));
b2_des=cross(b3_des,b1_des);
b2_des=b2_des/norm(b2_des);  %Normalized
R_des=[cross(b2_des,b3_des), b2_des, b3_des];
ang_err_mat=0.5*(R_des'*rot_mat-rot_mat'*R_des);
ang_err=[ang_err_mat(3,2), -ang_err_mat(3,1),ang_err_mat(2,1)]';
ang_rate_des=[0,0,yaw_rate_des(i)]';
ang_rate_err=ang_vel(:,i)-R_des'*rot_mat*ang_rate_des;
ang_moment=-att_gain*ang_err-ang_rate_gain*ang_rate_err+cross(ang_vel(:,i),I*ang_vel(:,i));
%ang_acc=-att_gain_mar*ang_err-ang_rate_gain_mar*ang_rate_err+cross(ang_vel(:,i),I*ang_vel(:,i));

%
thrust=-m*acc_W'*rot_mat(:,3);
%Convert to rotor velocities
omega_squared=inv_alloc_mat*[thrust;ang_moment];
%omega_squared=acc2omega*[thrust;ang_acc];
omega_squared(omega_squared<0)=0;
rotor_speeds(:,i)=sqrt(omega_squared);
%rotor_speeds(rotor_speeds(:,i)>omega_max,i)=omega_max;
if ~isreal(rotor_speeds(1,i))
    5
end
end
end