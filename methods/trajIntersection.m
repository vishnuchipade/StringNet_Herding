function [Ti]=trajIntersection(path1, path2, pathVel1, pathVel2, interSec)
% ------------------------------
% Author: Vishnu S. Chipade
% Email: vishnuc@umich.edu
% ------------------------------

%This function gives intersections segments between two trajectories as a
%time interval

global rho_D

%End points of the intersection segments
Si11=interSec.Si11;
Si12=interSec.Si12;
Si21=interSec.Si21;
Si22=interSec.Si22;

%Find the times at the end points of the intersection
Ti11=findTimeOnPath(Si11,path1,pathVel1);
Ti12=findTimeOnPath(Si12,path1,pathVel1);
Ti21=findTimeOnPath(Si21,path2,pathVel2);
Ti22=findTimeOnPath(Si22,path2,pathVel2);

if ~(Ti12<=Ti21 || Ti22<=Ti11)  %check if time intervals overlap
    %Find the overlapping time intervals
    Ti=[Ti11, Ti12, Ti21, Ti22];
    ascend_Ti=sort(Ti);
    %Pick the conflicting time interval
    Ti11_bar=ascend_Ti(2);
    Ti12_bar=ascend_Ti(3);
    %end points on path1 during this conflicting time interval
    Si11_bar=findPosOnPath(Ti11_bar,path1,pathVel1);
    Si12_bar=findPosOnPath(Ti12_bar,path1,pathVel1);
    %Pick some sample points on path1
    Si1_bar=Si11_bar:2*rho_D:Si12_bar;
    Si1_bar=[Si1_bar, Si12_bar];
    %Coordinates of the points in this intervals
    ri1_bar=findCoordOnPath(Si1_bar,path1);
    Ti1_bar=findTimeOnPath(Si1_bar,path1,pathVel1);
    
    nP=length(Si1_bar);
    
    Ti=[0,0]; %Zero norm vector, no collision
    i=0;
    while i<nP
        i=i+1;
        %Find the position of the robot 2 on path2 at T1_bar instances
        Si2_bar=findPosOnPath(Ti1_bar(i),path2,pathVel2);
        ri2_bar=findCoordOnPath(Si2_bar,path2);
        d=norm(ri1_bar(:,i)-ri2_bar);
        if d<=2*rho_D
            Ti(1)=Ti1_bar(i);  %first end point
            while d<=2*rho_D && i<nP
                i=i+1;
                Si2_bar=findPosOnPath(Ti1_bar(i),path2,pathVel2);
                ri2_bar=findCoordOnPath(Si2_bar,path2);
                d=norm(ri1_bar(:,i)-ri2_bar);
            end
            Ti(2)=Ti1_bar(i);  %second end point
        end
    end
else
    Ti=[0,0];  %No collision possible
end
end