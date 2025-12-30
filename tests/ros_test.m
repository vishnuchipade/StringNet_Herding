% ------------------------------
% Author: Vishnu S. Chipade
% Email: vishnuc@umich.edu
% ------------------------------

clear all
close all
rosshutdown

%% Master Initialisation
%rosinit;
%rosinit('http://dasc-2:11311/');
%rosinit('http://dasc-ThinkPad-P50:11311');
rosinit('http://192.168.1.119:11311');
%rosinit('http://localhost:11311')
%% Node initialisation
rovernumber = 4;
%R2_node = robotics.ros.Node('R2_node');
%R2_sub = robotics.ros.Subscriber(R2_node,'/vicon/R2/R2','geometry_msgs/TransformStamped');
%  formatspec = '/pelican/odometry_sensor1/pose';
%  string = sprintf(formatspec);
%  %RS_quad_pose = rossubscriber(string,'mav_msgs/Actuators');
%  R_sub_pose = rossubscriber(string,'geometry_msgs/Pose');
%RS_quad_pose = rossubscriber(string,'geometry_msgs/PoseStamped');

pause(5);
% pos = zeros(3,1);
% head = zeros(4,1);
% [time,pos,head] = statecallback(RS_quad_pose.LatestMessage,pos,head);

% formatspec2 = '/uav1/command/motor_speed';
% string2 = sprintf(formatspec2);
% R_pub_rot_speed = rospublisher(string2,'mav_msgs/Actuators');
% msg_pub_rot_speed = rosmessage(R_pub_rot_speed);
% msg_test.AngularVelocities= [590,600,600,650,600,600];
% send(R_pub,msg_test)
% R_sub_traj=rossubscriber('/firefly/command/trajectory','trajectory_msgs/MultiDOFJointTrajectory');
% %msg_traj=rosmessage(R_pub_traj);
% R_pub_traj=rospublisher('/firefly/command/trajectory','trajectory_msgs/MultiDOFJointTrajectory');
% msg_traj=rosmessage(R_pub_traj);
% R_pub_pose=rospublisher('/firefly/command/pose','geometry_msgs/PoseStamped');
% msg_pose=rosmessage(R_pub_pose);
%
% R_pub_traj=rospublisher('/firefly/command/roll_pitch_yawrate_thrust','mav_msgs/RollPitchYawrateThrust');
% msg_traj=rosmessage(R_pub_traj);
ND=5;
for i=1:ND
%Subcribe to odometry data
formatspec = '/defender%d/ground_truth/pose';
%formatspec = '/pelican/odometry_sensor1/pose';
string = sprintf(formatspec,i);
R_sub_pose(i) = rossubscriber(string,'geometry_msgs/Pose');

formatspec = '/defender%d/command/pose';
%formatspec = '/pelican/odometry_sensor1/pose';
string = sprintf(formatspec,i);
R_pub_pose (i)= rospublisher(string,'geometry_msgs/PoseStamped');
msg_pub_pose(i)=rosmessage(R_pub_pose(i));

formatspec = '/defender%d/ground_truth/odometry';
%formatspec = '/pelican/odometry_sensor1/odometry';
string = sprintf(formatspec,i);
R_sub_odometry(i) = rossubscriber(string,'nav_msgs/Odometry');

formatspec2 = '/defender%d/command/motor_speed';
string2 = sprintf(formatspec2,i);
R_pub_rot_speed(i)= rospublisher(string2,'mav_msgs/Actuators');
msg_pub_rot_speed(i) = rosmessage(R_pub_rot_speed(i));
end
client_gazebo_pause_physics = rossvcclient('/gazebo/pause_physics');
client_gazebo_unpause_physics = rossvcclient('/gazebo/unpause_physics');
rosTime=rossubscriber('/clock','rosgraph_msgs/Clock');
%%
omega=0.3;
R_circ=[4,4,4,4,6];
theta=[0,pi/2,pi,3*pi/2,pi/4];
rotor_speeds_history=[];
datum_time = rostime('now');
elapsedTime = 0;
dT=0;
counter = 0;
while(1)%elapsedTime<21)
    %msg_pose=receive(R_sub_pose);
    
    for i=1:ND
    
        pos_des=[R_circ(i)*cos(omega*dT+theta(i)),R_circ(i)*sin(omega*dT+theta(i)),.5*min(dT,30)]';
        vel_des=omega*[-R_circ(i)*sin(omega*dT+theta(i)),R_circ(i)*cos(omega*dT+theta(i)),.5*(sign(30-dT)+1)/2]';
        acc_des=-omega^2*[R_circ(i)*cos(omega*dT+theta(i)),R_circ(i)*sin(omega*dT+theta(i)),0]';
%     pos_des=[-10,1,1]';
%     vel_des=[0,0,0]';
%     acc_des=[0,0,0]';
    yaw_des=0;
    yaw_rate_des=0;
    %call(client_gazebo_unpause_physics);   %Unpause Gazebo physics engine
    
    msg_pose=R_sub_pose(i).LatestMessage;
    pos=[msg_pose.Position.X, msg_pose.Position.Y, msg_pose.Position.Z]';
    orientation=[msg_pose.Orientation.X, msg_pose.Orientation.Y, msg_pose.Orientation.Z,msg_pose.Orientation.W]';
    
    msg_odometry=R_sub_odometry(i).LatestMessage;
    Linear= msg_odometry.Twist.Twist.Linear;
    vel=[Linear.X,Linear.Y,Linear.Z]';
    Angular=msg_odometry.Twist.Twist.Angular;
    ang_vel=[Angular.X,Angular.Y,Angular.Z]';
    rotor_speeds=quad_controller(pos,vel,orientation,ang_vel,pos_des,vel_des,acc_des,yaw_des,yaw_rate_des);
    rotor_speeds_history{i}(:,counter+1)=rotor_speeds;
    
    %call(client_gazebo_unpause_physics);   %Unpause Gazebo physics engine
    %elapsedTime = rostime('now')-datum_time;
    %call(client_gazebo_pause_physics);  %Pause Gazebo physics engine
    elapsedTime = rosTime.LatestMessage.Clock_-datum_time; 
    %elapsedTime = rostime('now')-datum_time; 
    dT=double(elapsedTime.Sec)+double(elapsedTime.Nsec)*10^(-9);
    
    msg_pub_rot_speed(i).AngularVelocities=rotor_speeds;
    send(R_pub_rot_speed(i),msg_pub_rot_speed(i)) ;
    end
    %call(client_gazebo_pause_physics);  %Pause Gazebo physics engine
    counter = counter+1
    T(counter) = dT;
end
%msg_test.Angular.Z = 0;
%send(R2_pub,msg_test);
