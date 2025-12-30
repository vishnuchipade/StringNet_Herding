function [X1]=doubleIntDynamics(X0,U)
% ------------------------------
% Author: Vishnu S. Chipade
% Email: vishnuc@umich.edu
% ------------------------------

%Modified double integrator dynamics
%X=NA*[x y vx vy]'
%U=NA*[ux uy]'



N=length(X0)/4;
global dt;

[k1] = Xdot(X0,U);
[k2] = Xdot(X0 + dt/2*k1,U);
[k3] = Xdot(X0 + dt/2*k2,U);
[k4] = Xdot(X0 + dt*k3,U);

X1 = X0 + dt/6*(k1 + 2*k2 + 2*k3 + k4);


    function [f]=Xdot(X,U)
           
        f=zeros(4*N,1);
        for i=1:N
            f(4*i-3:4*i-2)=X(4*i-1:4*i);   %Only velocity is applied
            f(4*i-1:4*i)=U(2*i-1:2*i); 
        end
        
    end

end