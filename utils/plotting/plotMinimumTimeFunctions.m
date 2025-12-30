
%Minimum time for straight line segment under modified double integrator

% ------------------------------
% Author: Vishnu S. Chipade
% Email: vishnuc@umich.edu
% ------------------------------

%dynamics: s_dot=v; v_dot=u-Cd*v^2;


Gamma_s=.1;
vm=10;
Cd=0.1;
um=Cd*vm^2;

Np=200;
v1=0:vm/Np:vm;
v2=0:vm/Np:vm;

for i=1:Np+1
    
    for j=1:Np+1
        v0=v1(i);
        vf=v2(j);
        flag=0;
        if vf>v0
            G1=log((um-Cd*v0^2)/(um-Cd*vf^2))/(2*Cd);
            if Gamma_s>G1
                flag=1;
            end
        else
            G2=log((um+Cd*v0^2)/(um+Cd*vf^2))/(2*Cd);
            if Gamma_s>G2
                flag=1;
            end
        end
        
        if flag==1
        lambda=(um+Cd*vf^2)/(um-Cd*v0^2)*exp(2*Cd*Gamma_s);
        v_bar(i,j)=sqrt((lambda-1)/(lambda+1)*um/Cd);
        tau_s(i,j)=1/sqrt(um*Cd)*(2*atanh(sqrt(Cd*v_bar(i,j)^2/um))-atanh(sqrt(Cd*v0^2/um))-atanh(sqrt(Cd*vf^2/um)));
        else
            tau_s(i,j)=Inf;
        end
    end
end
figure
surface(v1,v2,tau_s);
xlabel('x');
ylabel('y');
zlabel('z')

%Find hessian of optimal t for v_dot=-Cd*v^2 + um-v^4/(R^2um); 
syms v0 vf G R vm um Cd
lambda=(um+Cd*vf^2)/(um-Cd*v0^2)*exp(2*Cd*G);
vbar=sqrt((lambda-1)/(lambda+1)*um/Cd);
t=sqrt(1/(Cd*um))*(2*atanh(sqrt(Cd*vbar^2/um))-atanh(sqrt(Cd*v0^2/um))-atanh(sqrt(Cd*vf^2/um)));
nabla2_t=hessian(t,[v0,vf]);
nabla2_t_fun=matlabFunction(nabla2_t, 'vars',[v0,vf,G,vm,um,Cd]);
G=40;
vm=10;
Cd=0.1;
um=Cd*vm^2;
vmc99=.99*vm;
v1=0:vmc99/100:vmc99;
v2=0:vmc99/100:vmc99;
 for i=1:length(v1)
    for j=1:length(v2)
        %vec=mat2cell([v1(i),v2(j),S,R,vm,um]);
        nabla2_t2{i,j}=nabla2_t_fun(v1(i),v2(j),G,vm,um,Cd);
        det_nabla2_t(i,j)=det(nabla2_t2{i,j});
    end
 end
 figure
 surf(v1,v2,det_nabla2_t)
 xlabel('x');
ylabel('y');
zlabel('z')
