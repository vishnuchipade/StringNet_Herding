%Exponential controller for DI dynamics 
% ------------------------------
% Author: Vishnu S. Chipade
% Email: vishnuc@umich.edu
% ------------------------------

%dynamics dot(r)=v; dot(v)=u
%controller u = sat_u1(-k*r) + sat_u2(-k*v)

%vector case i.e. r=[x,y], v=[vx,vy]
uMax=10;
uMax1=7;
uMax2=uMax-uMax1;
k1=1;

c1=0.12;
c2=5;
norm_r_bar = (uMax1/k1);
norm_v_bar= uMax2/k2;

global C_d dt
C_d=0.2;
dt=0.01;

A=[zeros(2),eye(2);-k1*eye(2),-k2*eye(2)];
Q=5*eye(4);
Ad=expm(A*dt);
P=lyap(A',Q)

V0_min=P(1,1)*norm_r_bar^2+P(3,3)*norm_v_bar^2-2*abs(P(1,3))*norm_r_bar*norm_v_bar
V0_max=P(1,1)*norm_r_bar^2+P(3,3)*norm_v_bar^2+2*abs(P(1,3))*norm_r_bar*norm_v_bar


NIter=50000;
NTheta=1;
for j=1:NTheta
    %theta=(j-1)*2*pi/(NTheta-1);
    theta=0;
    Theta(j)=theta;

X(:,1)=[405*rand(2,1);15*rand(2,1)];

time(1)=0;
traj_flag(1)=2;
for i=1:NIter
    norm_r=norm(X(1:2,i));
    norm_v=norm(X(3:4,i));
    if norm_r>1e-10
        phi = k1*X(1:2,i);
        norm_phi=norm(phi);
        U10=-(phi);
    else
        U10=[0,0]';
    end    
    if norm_v>1e-10 
        U20=-k2*(X(3:4,i));
    else
        U20=[0,0]';
    end
    norm_U10=norm(U10);
    norm_U20=norm(U20); 
    if norm_U10<1e-10
        U(:,i) =  min(uMax2,norm_U20)*U20/norm_U20;
    elseif norm_U20<1e-10
        U(:,i) = min(uMax1,norm_U10)*U10/norm_U10;
    else
    U(:,i) = min(uMax1,norm_U10)*U10/norm_U10 + min(uMax2,norm_U20)*U20/norm_U20;
    end
    
    X(:,i+1)=doubleIntDynamics(X(:,i),U(:,i));
    
    %Find the Lyapunov function at time i
    
     if norm_r<norm_r_bar && norm_v>norm_v_bar
        traj_flag(i)=1;
        %V(i)=1/alpha_prime*norm(phi)^alpha_prime+ c1* X(3:4,i)'*phi + c2*(uMax2*norm_v-V0_v); %*norm_v/uMax2
    elseif norm_r>norm_r_bar && norm_v>norm_v_bar
        traj_flag(i)=2;
        %V(i)=uMax1*norm_phi+ c1* X(3:4,i)'*phi + c2*(uMax2*norm_v-V0_v)-V0_phi;  %*norm_v/uMax2*norm_v/uMax1
    elseif norm_r>norm_r_bar && norm_v<norm_v_bar
        traj_flag(i)=3;
        %V(i)=uMax1*norm_phi+ c1* X(3:4,i)'*phi + c2/(3-alpha2)*norm_v^(3-alpha2)-V0_phi; %*norm_v/uMax1
    else
        traj_flag(i)=4;
        %V(i)=1/alpha_prime*norm(phi)^alpha_prime+ c1* X(3:4,i)'*phi + c2/(3-alpha2)*norm_v^(3-alpha2);
    end
    V(i)=X(:,i)'*P*X(:,i);
    
    f=[X(3:4,i);U(:,i)];
    %if i>1
        V_dot(i)=X(:,i)'*P*f+f'*P*X(:,i);
    %V_dot(i)=dV_by_dX'*f;
    %end
    V_dot_bar(i)=-X(:,i)'*Q*X(:,i);
    time(i+1)=time(i)+dt;
    if i>1
%     if V(i)-V(i-1)<1e-10 && V(i)<1e-15
%         break;
%     end
    end
end
V(i+1)=V(i);
V_dot(i+1)=V_dot(i);
V_dot_bar(i+1)=V_dot_bar(i);
convergence_time(j)=time(i)

end

colors={[1,0,0],[0,0,0],[0,0,1]};
markers={'--','d','*','.'};
circCos=cos(0:pi/100:2*pi);
circSin=sin(0:pi/100:2*pi);

figure('units','normalized','outerposition',[.1 0.1 .8 .8]);
R0=norm_r_bar;
plot(R0*circCos,R0*circSin)
hold on
for k=1:4
    ind=find(traj_flag==k);
plot(X(1,ind),X(2,ind),markers{k},'color',colors{1},'markersize',3)
hold on
end

figure('units','normalized','outerposition',[.1 0.1 .8 .8]);
for k=1:4
    ind=find(traj_flag==k);
plot(time(ind),X(1,ind),markers{k},'color',colors{1},'markersize',3)
hold on
plot(time(ind),X(2,ind),markers{k},'color',colors{2},'markersize',3)
end


figure('units','normalized','outerposition',[.1 0.1 .8 .8]);
for k=1:4
    ind=find(traj_flag==k);
    plot(time(ind),V_dot(ind),markers{k},'color',colors{1},'markersize',3)
hold on
plot(time(ind),V_dot_bar(ind),markers{k},'color',colors{2},'markersize',3)
 plot(time(ind),V(ind),markers{k},'color',colors{3},'markersize',3)
end
plot(time,V0_min*ones(size(time)),'g')
plot(time,V0_max*ones(size(time)),'g')
xlabel('t')
ylabel('V')
legend('$\dot{V}$','$\bar{\dot{V}}$','-V')


figure
%subplot(2,1,1)
plot(time,V)
xlabel('t')
ylabel('V')

figure
plot(Theta,convergence_time)
xlabel('t')
ylabel('Convergence Time')