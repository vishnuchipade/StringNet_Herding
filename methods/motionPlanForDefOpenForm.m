function motionPlan=motionPlanForDefOpenForm(tanG,XD,XD_des,ND,flagPlotPaths,figNumber)
% ------------------------------
% Author: Vishnu S. Chipade
% Email: vishnuc@umich.edu
% ------------------------------

global v_maxDC u_maxD v_maxD C_d
%Find the shortest paths and speed profiles

umd=0.8*u_maxD(1);
vmd=sqrt(umd/C_d);
for jj=1:ND
    for j=1:ND
        %Column represents the defender and the row goal location
        [Path{jj,j},tanG_prime{jj,j}]=findShortestPath(tanG,XD(1:2,j),XD_des(1:2,jj));
        if (flagPlotPaths)
        PlotShortestPath(tanG,Path{jj,j},tanG_prime{jj,j},figNumber,[.5,0.5,0.5],'--',1.2)
        end
        pathVel{jj,j}=findPathSpeeds(Path{jj,j},v_maxDC(j), umd, vmd);  %rI=XD(1:2,j), rF=XD_des(1:2,jj)
        optT(jj,j)=pathVel{jj,j}.T(end);
    end
end


Path0=Path(:);  %arrange the paths such paths of D1 are the first and then that of D2 and so on
%Plot the numbers corrsponding to the defenders
if(0)
for j=1:ND
    text(XD(1,j),XD(2,j),['$\mathcal{D}_', num2str(j), '$'])
end
end
%%
Pbar=diag(zeros(1,ND^2));
for j=1:ND^2
    for jj=j+1:ND^2
        interSec{j,jj}=pathIntersections(Path0{j},Path0{jj});
        Pbar(j,jj)=interSec{j,jj}.Pbar1;
        Pbar(jj,j)=interSec{j,jj}.Pbar2;
    end
end
tic
[assign,assignCost,maxOptT,results,model]=defGoalAssignMIQP(optT(:),Pbar,ND);
tCompMIQP=toc;
%%
%     assign=[2,1,3];
%Assigned path and related values
for j=1:ND
    assignedPath{j}=Path{assign(j),j};
    assignedTanG_prime{j}=tanG_prime{assign(j),j};
    assignedPathVel{j}=pathVel{assign(j),j};
    assignedOptT(j)=pathVel{assign(j),j}.T(end);
    if (0)
    PlotShortestPath(tanG,Path{assign(j),j},[0.1,0,0.75],tanG_prime{assign(j),j},2);
    end
end

for j=1:ND
    jp=(j-1)*ND+assign(j);
    for jj=j+1:ND
        jjp=(jj-1)*ND+assign(jj);
        assignedInterSec{j,jj}=interSec{jp,jjp};
    end
end
%%
%Find the terminal times at the intersection segments on the
%assigned paths

inputStartTime.optT=assignedOptT;
NiSeg=0;  %Number of intersection segments on all the paths
leadTime=[];
for j=1:ND
    for jj=j+1:ND
        if assignedInterSec{j,jj}.flag0==1
            NiSeg=NiSeg+1;
            S11=assignedInterSec{j,jj}.Si11;
            S12=assignedInterSec{j,jj}.Si12;
            assignedInterSec{j,jj}.T11=findTimeOnPath(S11,assignedPath{j},assignedPathVel{j});
            assignedInterSec{j,jj}.T12=findTimeOnPath(S12,assignedPath{j},assignedPathVel{j});
            S21=assignedInterSec{j,jj}.Si21;
            S22=assignedInterSec{j,jj}.Si22;
            assignedInterSec{j,jj}.T21=findTimeOnPath(S21,assignedPath{jj},assignedPathVel{jj});
            assignedInterSec{j,jj}.T22=findTimeOnPath(S22,assignedPath{jj},assignedPathVel{jj});
            %interSecTime(NiSeg).T11=assignedInterSec{j,jj}.T11;
            leadTime{j,jj}.T1=leadTimeToAvoidCollision(assignedPath{j},assignedPath{jj},assignedPathVel{j},assignedPathVel{jj},assignedInterSec{j,jj},1);
            leadTime{j,jj}.T2=leadTimeToAvoidCollision(assignedPath{j},assignedPath{jj},assignedPathVel{j},assignedPathVel{jj},assignedInterSec{j,jj},2);
        end
    end
end
inputStartTime.NiSeg=NiSeg;
inputStartTime.interSec=assignedInterSec;
inputStartTime.leadTime=leadTime;

%Find the optimal start time schedule for all the defenders
tic
if ~isempty(leadTime)
startTime=startTimeForCollisionAvoidance(inputStartTime);
else
 startTime=zeros(ND,1);
end
tCompMILP=toc; %Time to calculate the lead times and solve the MILP

%club together the required results and output them
motionPlan.assign=assign;
motionPlan.assignCost=assignCost;
motionPlan.maxOptT=maxOptT;
motionPlan.Path=assignedPath;
motionPlan.tanG_prime=assignedTanG_prime;
motionPlan.interSec=assignedInterSec;
motionPlan.pathVel=assignedPathVel;
motionPlan.optT=assignedOptT;
motionPlan.leadTime=leadTime;
motionPlan.startTime=startTime;
motionPlan.tCompMIQP=tCompMIQP;
motionPlan.tCompMILP=tCompMILP;
motionPlan.tComp=tCompMIQP+tCompMILP;
motionPlan.XD=XD;
motionPlan.XD_des=XD_des;
end