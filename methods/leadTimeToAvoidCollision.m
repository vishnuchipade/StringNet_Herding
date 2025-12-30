function leadTime=leadTimeToAvoidCollision(path1,path2,pathVel1,pathVel2,interSec,leadId)
% ------------------------------
% Author: Vishnu S. Chipade
% Email: vishnuc@umich.edu
% ------------------------------

%This function calculates the lead time on path1 to avoid collision with
%the robot2 on path2 with given velocity profiles and the intersection
%segments

%The 'lead time' defines the time 'robot1' should lead 'robot2' at the beginning
%of the intersection segment 'interSec'

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

if leadId==1   %Robot 1 is leading 
    
    %Check if the time intervals corresponding to the collision segments overlap
    if Ti12>Ti21
        pathVel20=pathVel2;
        %intervals of time for bisection method;
        t1=0;
        t2=Ti12-Ti21;
        Ti0=trajIntersection(path1,path2,pathVel1,pathVel20,interSec);
        if ~isempty(Ti0)
            leadTime=t2-t1;
            pathVel20=pathVel2;
            pathVel20.T=pathVel20.T+leadTime;
            pathVel20.T_bar1=pathVel20.T_bar1+leadTime;
            pathVel20.T_bar2=pathVel20.T_bar2+leadTime;
            Ti=trajIntersection(path1,path2,pathVel1,pathVel20,interSec);
            while (norm(Ti-Ti0)>1e-5 || norm(Ti)<1e-5) && leadTime>1e-5
                dTi0=Ti0(2)-Ti0(1);
                dTi=Ti(2)-Ti(1);
                if dTi<=dTi0
                    t2=(t1+t2)/2;
                    Ti0=Ti;
                    
                    leadTime=t2-t1;
                    pathVel20=pathVel2;
                    pathVel20.T=pathVel20.T+leadTime;
                    pathVel20.T_bar1=pathVel20.T_bar1+leadTime;
                    pathVel20.T_bar2=pathVel20.T_bar2+leadTime;
                    Ti=trajIntersection(path1,path2,pathVel1,pathVel20,interSec);
                elseif dTi>dTi0
                    t1=(t1+t2)/2;
                    Ti0=Ti;
                    
                    leadTime=t2-t1;
                    pathVel20=pathVel2;
                    pathVel20.T=pathVel20.T+leadTime;
                    pathVel20.T_bar1=pathVel20.T_bar1+leadTime;
                    pathVel20.T_bar2=pathVel20.T_bar2+leadTime;
                    Ti=trajIntersection(path1,path2,pathVel1,pathVel20,interSec);
                end
            end
        else
            leadTime=0;
        end
    else
        leadTime=0;
    end
    
else %leadId==2  when robot 2 is leading
    
    %Check if the time intervals corresponding to the collision segments overlap
    if Ti22>Ti11
        pathVel10=pathVel1;
        %intervals of time for bisection method;
        t1=0;
        t2=Ti22-Ti11;
        Ti0=trajIntersection(path1,path2,pathVel10,pathVel2,interSec);
        if ~isempty(Ti0)
            leadTime=t2-t1;
            pathVel10=pathVel1;
            pathVel10.T=pathVel10.T+leadTime;
            pathVel10.T_bar1=pathVel10.T_bar1+leadTime;
            pathVel10.T_bar2=pathVel10.T_bar2+leadTime;
            Ti=trajIntersection(path1,path2,pathVel10,pathVel2,interSec);
            while (norm(Ti-Ti0)>1e-5 || norm(Ti)<1e-5) && leadTime>1e-5
                dTi0=Ti0(2)-Ti0(1);
                dTi=Ti(2)-Ti(1);
                if dTi<=dTi0
                    t2=(t1+t2)/2;
                    Ti0=Ti;
                    
                    leadTime=t2-t1;
                    pathVel10=pathVel1;
                    pathVel10.T=pathVel10.T+leadTime;
                    pathVel10.T_bar1=pathVel10.T_bar1+leadTime;
                    pathVel10.T_bar2=pathVel10.T_bar2+leadTime;
                    Ti=trajIntersection(path1,path2,pathVel10,pathVel2,interSec);
                elseif dTi>dTi0
                    t1=(t1+t2)/2;
                    Ti0=Ti;
                    
                    leadTime=t2-t1;
                    pathVel10=pathVel1;
                    pathVel10.T=pathVel10.T+leadTime;
                    pathVel10.T_bar1=pathVel10.T_bar1+leadTime;
                    pathVel10.T_bar2=pathVel10.T_bar2+leadTime;
                    Ti=trajIntersection(path1,path2,pathVel10,pathVel2,interSec);
                end
            end
        else
            leadTime=0;
        end
    else
        leadTime=0;
    end
end

end