function [lambda]=lambdaInterSecCircLine(rc,rhoc,r1,mL,cL,drx,dry)
% ------------------------------
% Author: Vishnu S. Chipade
% Email: vishnuc@umich.edu
% ------------------------------

global rho_D
if mL<1e16
    b=-2*(rc(1)-mL*cL+mL*rc(2));
    a=(1+mL^2);
    c=rc(1)^2+cL^2+rc(2)^2-2*cL*rc(2)-(rhoc)^2;
    x_int01=(-b-sqrt(b^2-4*a*c))/(2*a);
    %y_int01=mL1*x_int+cL1;
    lambda01=(x_int01-r1(1))/drx;
    x_int02=(-b+sqrt(b^2-4*a*c))/(2*a);
    %y_int02=mL1*x_int02+cL1;
    lambda02=(x_int02-r1(1))/drx;
else
    x_int01=r1(1);
    y_int01=rc(2)+sqrt((rhoc)^2-(x_int01-rc(1))^2);
    y_int02=rc(2)-sqrt((rhoc)^2-(x_int01-rc(1))^2);
    lambda01=(y_int01-r1(2))/dry;
    lambda02=(y_int02-r1(2))/dry;
end
lambda=[lambda01,lambda02];
% x_int=[x_int1,x_int2];
% y_int=[y_int1,y_int2];
end