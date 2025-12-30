function [startTime] = startTimeForCollisionAvoidance(inpStruct)
% ------------------------------
% Author: Vishnu S. Chipade
% Email: vishnuc@umich.edu
% ------------------------------

%This function gives initial times for each robot to start executing their
%motion on given paths with given velocity profiles such that no two robot
%collide, this is based on sufficient conditions only i.e. no two robots
%are inside their collision segments at the same. This could be
%conservative in certain cases where there is a significant overlap between
%two paths. The improved version with necessary conditions is developed in
%other function

%Input: 'inpStruct', struct with intersection segments and corresponding time
%interval information
%Output: startTime, an array with optimal starting times


%MILP formulation

%Number of decision variabls
%x=[t_m, t_s1,t_s2,...t_sND, delta_1, delta_2,..., delta_NiSeg]
%Nx=1+ND+2*NiSeg
NiSeg=inpStruct.NiSeg;
interSec=inpStruct.interSec;
optT=inpStruct.optT;

ND=size(interSec,2);
Nx=1+ND+NiSeg;  %number of decision variables
Nc=ND+2*NiSeg+ND; %number of constraints


%optT=interSec.optT;
M=sum(optT); %big number for linear constraints using binary variables
%c vector
model.obj=zeros(Nx,1);
model.obj(1)=1;

%Constraints
A=zeros(Nc,Nx);
b=zeros(Nc,1);

%objective constraints for representing the problem as LP
for j=1:ND
    A(j,1)=-1;
    A(j,1+j)=1;
    b(j,1)=-optT(j);
end

%collision constraints
ci0=1+ND;
ci=0; %index of the current binary variabel under consideration (i.e. current intersection segment)
for j=1:ND
    for jj=j+1:ND
        if interSec{j,jj}.flag0==1
            ci=ci+1;
            A(ND+2*ci-1, ci0+ci)=-M;
            A(ND+2*ci-1, 1+j)=1;
            A(ND+2*ci-1, 1+jj)=-1;
            b(ND+2*ci-1,1)= interSec{j,jj}.T21-interSec{j,jj}.T12;
            A(ND+2*ci, ci0+ci)=M;
            A(ND+2*ci, 1+j)=-1;
            A(ND+2*ci, 1+jj)=1;
            b(ND+2*ci,1)= -interSec{j,jj}.T22+interSec{j,jj}.T11+M;            
        end
    end    
end

ci02=ND+2*ci;
%feasibility constraints
for j=1:ND
    A(ci02+j,1+j)=-1;
end


model.A=sparse(A);
model.rhs=b;
model.sense = '<';
model.vtype(1:ND+1) = 'C';
model.vtype(ND+2:Nx) = 'B';

gurobi_write(model, 'startTimeForCollisionAvoidance.lp');

results = gurobi(model);
startTime=results.x(2:ND+1);


end

