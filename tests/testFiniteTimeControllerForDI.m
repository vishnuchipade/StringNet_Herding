%Fintie time controller for DI dynamics 
% ------------------------------
% Author: Vishnu S. Chipade
% Email: vishnuc@umich.edu
% ------------------------------

%dynamics dot(r)=v; dot(v)=u
%controller u = sat_u1(-k2*phi*||phi||^(alpha1-1)) + sat_u2(-k2*v*||v||^(alpha2-1)
%where phi= r + 1/(k2(2-alpha2))*v*||v||^(1-alpha2)

%scalar case i.e. r=[x], v=[vx]
uMax=1000;
uMax1=700;
uMax2=uMax-uMax1;
k1=1;
k2=1;
c1=0.2;
c2=3;
alpha2=0.4;
alpha1=alpha2/(2-alpha2);
alpha_prime=(3-alpha2)/(2-alpha2);
norm_r_bar = (uMax1/k1)^(1/alpha1);
V0=k1*alpha1/(1+alpha1)*(uMax1/k1)^((1+alpha1)/alpha1);
global C_d dt
C_d=0.2;
dt=0.0001;
NIter=50000;
NTheta=2;
for j=1:NTheta
    theta=(j-1)*2*pi/(NTheta-1);
    Theta(j)=theta;

X(:,1)=rand(2,1);

time(1)=0;
for i=1:NIter
    norm_r=norm(X(1,i));
    norm_v=norm(X(2,i));
    if norm_r>1e-10
        phi = X(1,i) + alpha1/alpha2*X(2,i)*norm_v^(1-alpha2);
        U10=-(phi)*norm(phi)^(alpha1-1);
    else
        U10=[0,0]';
    end    
    if norm_v>1e-10 
        U20=-(X(2,i))*norm_v^(alpha2-1);
    else
        U20=[0]';
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
    
    X(:,i+1)=X(:,i)+dt*[X(2,i),U(:,i)]';
    
    %Find the Lyapunov function at time i
    
    V(i)=1/alpha_prime*norm(phi)^alpha_prime+ c1* X(2,i)'*phi + c2/(3-alpha2)*norm_v^(3-alpha2);
    
    f=[X(2,i);U(:,i)];
    if i>1
        V_dot(i)=(V(i)-V(i-1))/dt;
    %V_dot(i)=dV_by_dX'*f;
    end
    V_dot_bar(i)=-0.1*V(i)^(2/(3-alpha2));
    time(i+1)=time(i)+dt;
    if i>1
    if V(i)-V(i-1)<1e-10 && V(i)<1e-5
        break;
    end
    end
end
V(i+1)=V(i);
V_dot(i+1)=V_dot(i);
V_dot_bar(i+1)=V_dot_bar(i);
convergence_time(j)=time(i)

end

figure
plot(time,X(1,:),time,X(2,:))

figure
%subplot(2,1,1)
plot(time,V)
xlabel('t')
ylabel('V')

figure
plot(time,V_dot,time,V_dot_bar,time,-V)%,time,V1_dot,time,V2_dot)
xlabel('t')
ylabel('V')
legend('$\dot{V}$','$\bar{\dot{V}}$','-V')

figure
plot(time,V_dot_bar)


figure
plot(Theta,convergence_time)
xlabel('t')
ylabel('Convergence Time')



%%

%Fintie time controller for DI dynamics 
%dynamics dot(r)=v; dot(v)=u
%controller u = sat_u1(-k2*phi*||phi||^(alpha1-1)) + sat_u2(-k2*v*||v||^(alpha2-1)
%where phi= r + 1/(k2(2-alpha2))*v*||v||^(1-alpha2)

%vector case i.e. r=[x,y], v=[vx,vy]
uMax=10;
uMax1=7;
uMax2=uMax-uMax1;
k1=1;
k2=1;
c1=0.12;
c2=5;
alpha2=0.7;
alpha1=alpha2/(2-alpha2);
alpha_prime=(3-alpha2)/(2-alpha2);
norm_r_bar = (uMax1/k1)^(1/alpha1);
norm_phi_bar = (uMax1)^(1/alpha1);
norm_v_bar=uMax2^(1/alpha2);
V0=k1*alpha1/(1+alpha1)*(uMax1/k1)^((1+alpha1)/alpha1);
V0_phi=uMax1*norm_phi_bar-1/alpha_prime*norm_phi_bar^alpha_prime;
V0_v=uMax2*norm_v_bar-1/(3-alpha2)*norm_v_bar^(3-alpha2);
global C_d dt
C_d=0.2;
dt=0.01;
NIter=50000;
NTheta=2;
for j=1:NTheta
    theta=(j-1)*2*pi/(NTheta-1);
    Theta(j)=theta;

X(:,1)=[405*rand(2,1);15*rand(2,1)];

time(1)=0;
traj_flag(1)=2;
for i=1:NIter
    norm_r=norm(X(1:2,i));
    norm_v=norm(X(3:4,i));
    if norm_r>1e-10
        phi = X(1:2,i) + alpha1/alpha2*X(3:4,i)*norm_v^(1-alpha2);
        norm_phi=norm(phi);
        U10=-(phi)*norm(phi)^(alpha1-1);
    else
        U10=[0,0]';
    end    
    if norm_v>1e-10 
        U20=-(X(3:4,i))*norm_v^(alpha2-1);
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
    
    V(i)=1/alpha_prime*norm(phi)^alpha_prime+ c1* X(3:4,i)'*phi + c2/(3-alpha2)*norm_v^(3-alpha2);
    
    if norm_phi<norm_phi_bar && norm_v>norm_v_bar
        traj_flag(i)=1;
        V(i)=1/alpha_prime*norm(phi)^alpha_prime+ c1* X(3:4,i)'*phi + c2*(uMax2*norm_v-V0_v); %*norm_v/uMax2
    elseif norm_phi>norm_phi_bar && norm_v>norm_v_bar
        traj_flag(i)=2;
        V(i)=uMax1*norm_phi+ c1* X(3:4,i)'*phi + c2*(uMax2*norm_v-V0_v)-V0_phi;  %*norm_v/uMax2*norm_v/uMax1
    elseif norm_phi>norm_phi_bar && norm_v<norm_v_bar
        traj_flag(i)=3;
        V(i)=uMax1*norm_phi+ c1* X(3:4,i)'*phi + c2/(3-alpha2)*norm_v^(3-alpha2)-V0_phi; %*norm_v/uMax1
    else
        traj_flag(i)=4;
        V(i)=1/alpha_prime*norm(phi)^alpha_prime+ c1* X(3:4,i)'*phi + c2/(3-alpha2)*norm_v^(3-alpha2);
    end
    
    
    f=[X(3:4,i);U(:,i)];
    if i>1
        V_dot(i)=(V(i)-V(i-1))/dt;
    %V_dot(i)=dV_by_dX'*f;
    end
    V_dot_bar(i)=-0.2*V(i)^(2/(3-alpha2));
    time(i+1)=time(i)+dt;
    if i>1
    if V(i)-V(i-1)<1e-10 && V(i)<1e-15
        break;
    end
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
R0=norm_phi_bar-1/(2-alpha2)*uMax^((2-alpha2)/alpha2);
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
