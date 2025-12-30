function [model,results]=test_qp()
% ------------------------------
% Author: Vishnu S. Chipade
% Email: vishnuc@umich.edu
% ------------------------------

% Copyright 2019, Gurobi Optimization, LLC
%
% This example formulates and solves the following simple QP model:
%  minimize
%      x^2 + x*y + y^2 + y*z + z^2 + 2 x
%  subject to
%      x + 2 y + 3 z >= 4
%      x +   y       >= 1
%      x, y, z non-negative
%
% It solves it once as a continuous model, and once as an integer
% model.

%names = {'x1', 'x2', 'x3', 'x4', 'x5','x6', 'x7', 'x8', 'x9'};
names = {'x1', 'x2', 'x3', 'x4'};
model.varnames = names;
% Q=zeros(9,9);
% Q(1,4:9)=[1,0,0,0,2,0];
% Q(2,4:9)=[5,2,0,1,0,0];
% Q(3,4:9)=[0,2,0,3,0,1];
% Q(4,[1:3,7:9])=[0,2,0,3,0,1];
% Q(5,[1:3,7:9])=[5,2,0,1,0,0];
% Q(6,[1:3,7:9])=[1,0,0,0,2,0];
% Q(7,1:6)=[5,2,0,1,0,0];
% Q(8,1:6)=[1,2,0,3,0,1];
% Q(9,1:6)=[1,0,0,0,2,0];
% 
% model.Q = sparse(Q);
% c=rand(9,1);
% model.c=c;
% 
% A=zeros(3,9);
% A(1,1:3)=ones(1,3);
% A(2,4:6)=ones(1,3);
% A(3,7:9)=ones(1,3);
% A(4,1:3:9)=ones(1,3);
% A(5,2:3:9)=ones(1,3);
% A(6,3:3:9)=ones(1,3);
% b=ones(6,1);
% model.A = sparse(A);
% model.b=b;

Q=diag(ones(4,1));
Q(2,3)=1;
Q(3,2)=1;
Q(1,4)=1;
Q(4,1)=1;

model.Q = sparse(Q);
c=rand(4,1);
model.c=c;

A=zeros(4,4);
A(1,1:2)=ones(1,2);
A(2,3:4)=ones(1,2);
A(3,[1,3])=ones(1,2);
A(4,[2,4])=ones(1,2);
b=ones(4,1);
model.A = sparse(A);
%model.b=b;


model.obj = c;
model.rhs = b;
model.sense = '=';

gurobi_write(model, 'qp.lp');

results = gurobi(model);

for v=1:length(names)
    fprintf('%s %e\n', names{v}, results.x(v));
end

fprintf('Obj: %e\n', results.objval);

model.vtype = 'B';
model.isQP=1;
results  = gurobi(model);

for v=1:length(names)
    fprintf('%s %e\n', names{v}, results.x(v));
end

fprintf('Obj: %e\n', results.objval);

end
