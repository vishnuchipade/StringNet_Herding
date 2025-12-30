function [assign,cost,maxOptT, results,model]=defGoalAssignMIQP(optT,Pbar,ND)
% ------------------------------
% Author: Vishnu S. Chipade
% Email: vishnuc@umich.edu
% ------------------------------


%objective function
model.Q = sparse(Pbar); %quadratic matrix
model.obj=optT; % linear vector

%Constraints
for j=1:ND
    A(j,ND*(j-1)+1:ND*j)=ones(1,ND);
    A(ND+j,[j:ND:ND^2])=ones(1,ND);
end
model.A=sparse(A);
b=ones(2*ND,1);
model.rhs=b;
model.sense = '=';
model.vtype = 'B';

gurobi_write(model, 'defGoalAssignMIQP.lp');

results = gurobi(model);

%optimal cost
cost=model.obj'*results.x+results.x'*model.Q*results.x;
%maximum optimal time to reach all the defenders
maxOptT=max(model.obj.*results.x);

%find the assignment from the solution
tempInd=find(results.x==1);
assign=mod(tempInd,ND);
assign(assign==0)=ND;
end