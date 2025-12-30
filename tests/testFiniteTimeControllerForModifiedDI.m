%Fintie time controller for DI dynamics with quadratic drag term
% ------------------------------
% Author: Vishnu S. Chipade
% Email: vishnuc@umich.edu
% ------------------------------

%controller u = sat_u1(-k1*phi*||phi||^(alpha1-1)) + sat_u2(-k2*v*||v||^(alpha2-1)+C_d*||v||v)
%where phi= r + 1/(k2(2-alpha2))*v*||v||^(1-alpha2)
%first slows down the agent to zero velocity and then starts the whole
%controller (we can find the bound on time of convergence analytically)
%The nominal finite time controller withour explicit gaurantees on time of
%convregence is simualted later in the script


uMax=10;
uMax1=7;
uMax2=uMax-uMax1;


alpha2=0.9;
alpha1=alpha2/(2-alpha2);
alpha_prime=(3-alpha2)/(2-alpha2);

global C_d dt
C_d=0.2;

v_max=sqrt(uMax/C_d);

k2=(uMax-uMax2)/v_max^(alpha2);
k1=k2;

f=@(x) C_d*x^2-k2*x^alpha2-uMax2;
norm_v_bar = abs(fsolve(f,sqrt(uMax2/C_d)));
norm_v_bar_alpha=(uMax1/k1)^(1/alpha2);
m1=-k2*alpha2*norm_v_bar_alpha^(alpha2-1);
f=@(x) (alpha2^2-1)*k2^2*x^(2*alpha2)+2*uMax1*k2*x^(alpha2)-uMax1^2;
% v_tilde=fsolve(f,norm_v_bar_alpha);
% m1_tilde=-k2*alpha2*v_tilde^(alpha2-1);
% a_tilde=uMax1-k2*v_tilde^alpha2;
% t_v_tilde=1/m1_tilde*log(a_tilde/(a_tilde-m1_tilde*v_tilde));
% Delta=(a_tilde-m1_tilde*v_tilde)/m1_tilde;
% r_v_tilde=Delta*exp(m1_tilde*t_v_tilde)/m1_tilde-Delta*t_v_tilde-Delta;
N_vl0=1000;

vl0=norm_v_bar_alpha/N_vl0:norm_v_bar_alpha/N_vl0:norm_v_bar_alpha;

if(1)
for j=1:length(vl0)-1    
    if j==1
        m_vl0(j)=-k2*alpha2*vl0(j)^(alpha2-1);
        m_vl0(j+1)=-k2*alpha2*vl0(j+1)^(alpha2-1);
        al0(j)=uMax1-k2*vl0(j)^alpha2;
        al0(j+1)=uMax1-k2*vl0(j+1)^alpha2;
        c0(j)=al0(j)-m_vl0(j)*vl0(j);
        c0(j+1)=al0(j+1)-m_vl0(j+1)*vl0(j+1);
        vli(j+1)= (c0(j+1)-c0(j))/(m_vl0(j)-m_vl0(j+1));
        ali(j+1)=m_vl0(j)*vli(j+1)+c0(j);
        tli(j+1)=log((al0(j)-m_vl0(j)*vl0(j)+m_vl0(j)*vli(j+1))/(al0(j)-m_vl0(j)*vl0(j)))/m_vl0(j);
        Delta(j)=(al0(j)-m_vl0(j)*vl0(j))/m_vl0(j);
        rli(j+1)=Delta(j)/m_vl0(j)*exp(m_vl0(j)*tli(j+1))-Delta(j)*tli(j+1)-Delta(j)/m_vl0(j);
    else
        m_vl0(j+1)=-k2*alpha2*vl0(j+1)^(alpha2-1);
        al0(j+1)=uMax1-k2*vl0(j+1)^alpha2;
        c0(j+1)=al0(j+1)-m_vl0(j+1)*vl0(j+1);
        vli(j+1)= (c0(j+1)-c0(j))/(m_vl0(j)-m_vl0(j+1));
        ali(j+1)=m_vl0(j)*vli(j+1)+c0(j);
        tli(j+1)=log((al0(j)-m_vl0(j)*vl0(j)+m_vl0(j)*vli(j+1))/(al0(j)-m_vl0(j)*vl0(j)+m_vl0(j)*vli(j)))/m_vl0(j);
        Delta(j)=(al0(j)-m_vl0(j)*vl0(j))/m_vl0(j);
        rli(j+1)=rli(j)+(Delta(j)+vli(j))/m_vl0(j)*exp(m_vl0(j)*tli(j+1))-Delta(j)*tli(j+1)-(Delta(j)+vli(j))/m_vl0(j);
    end
end
tli(1)=0;
vli(1)=0;
ali(1)=c0(1);
rli(1)=0;
j=j+1;
Delta(j)=(al0(j)-m_vl0(j)*vl0(j))/m_vl0(j);
end

figure
v=0:norm_v_bar_alpha/100:norm_v_bar_alpha;
plot(v,uMax1-k2*v.^alpha2,'b',vli,ali,'r')


norm_r_bar = ((uMax1/k1)^(1/alpha1)-1/(k2*(2-alpha2))*v_max^(2-alpha2));
norm_phi_bar = (uMax1/k1)^(1/alpha1);

norm_r_bar2 = ((uMax1/k1)^(1/alpha1)+1/(k2*(2-alpha2))*v_max^(2-alpha2));


V0=k1*alpha1/(1+alpha1)*(uMax1/k1)^((1+alpha1)/alpha1);
V0_phi=uMax1*norm_phi_bar-1/alpha_prime*norm_phi_bar^alpha_prime;
V0_v=uMax2*norm_v_bar-1/(3-alpha2)*norm_v_bar^(3-alpha2);

c1=0.99;
c2=1.1;

figure
v=0:norm_v_bar_alpha/100:norm_v_bar_alpha;
plot(v,c2*(uMax1.*v.^(2-alpha2)-k2*v.^2)-uMax1^((5-2*alpha2)/(2-alpha2))/k2*v.^(1-alpha2)+c1*uMax1/k2*v.^(2-alpha2)-c1*uMax1+c1*k2*uMax1*v.^alpha2)


c_bar=0.5;
dt=0.01;
NIter=50000;
NTheta=1;
speed0=v_max;
t_norm_r_bar=0;
t_norm_phi_bar=0;
for j=1:NTheta
     %theta=(j-1)*2*pi/(NTheta-1);
     theta=pi/3;
     Theta(j)=theta;

X(:,1)=[norm_r_bar+0.3;0;speed0*cos(theta);speed0*sin(theta)];
%X(:,1)=[40*rand(2,1);-7*rand(2,1)];
X0=X(:,1);
norm_r0=norm(X(1:2,1));
norm_v0=norm(X(3:4,1));

flag_v=0;
time(1)=0;
for i=1:NIter
    norm_r=norm(X(1:2,i));
    norm_v=norm(X(3:4,i));
    if norm_r>1e-10
        phi =X(1:2,i)+ 1/(k2*(2-alpha2))*X(3:4,i)*norm_v^(1-alpha2);
        norm_phi=norm(phi);
        U10=-k1*(phi)*norm(phi)^(alpha1-1);
    else
        U10=[0,0]';
    end
    
    norm_v_arr(i)=norm_v;
    if norm_v>1e-10 
        U20=-k2*(X(3:4,i))*norm(X(3:4,i))^(alpha2-1)+C_d*norm(X(3:4,i))*X(3:4,i);
    else
        U20=[0,0]';
    end
    norm_U10=norm(U10);
    norm_U20=norm(U20);
    if norm_v>1e-5 && flag_v==0
    t_v0=time(i);
    else
    flag_v=1;
    end
    uMax1_prime=uMax1;%*(min(1,time(i)-t_v0));
    if norm_U10<1e-10
        U(:,i) =  min(uMax2,norm_U20)*U20/norm_U20;
    elseif norm_U20<1e-10
        U(:,i) = min(uMax1_prime,norm_U10)*U10/norm_U10;
    else
    U(:,i) = min(uMax1_prime,norm_U10)*U10/norm_U10 + min(uMax2,norm_U20)*U20/norm_U20;
    end
    norm_U_arr(i)=norm(U(:,i));
    X(:,i+1)=modifiedDIDynamics(X(:,i),U(:,i));
    
    %Find the Lyapunov function at time i
    
    if norm_r<norm_r_bar
        V(i)=k1/(1+alpha1)*norm_r^(1+alpha1) + 0.5*norm_v^2;
        dV_by_dX=[k1*X(1:2,i)*norm_r^(alpha1-1);X(3:4,i)];
    else
        V(i)=uMax1*norm_r + 0.5*norm_v^2 -V0;
        dV_by_dX=[uMax*X(1:2,i)*norm_r;X(3:4,i)];
    end
    
     if norm_phi<norm_phi_bar && norm_v>norm_v_bar
        traj_flag(i)=1;
        V(i)=1/alpha_prime*norm(phi)^alpha_prime+ c1* X(3:4,i)'*phi + c2*(uMax2*norm_v-V0_v); %*norm_v/uMax2
    elseif norm_phi>norm_phi_bar && norm_v>norm_v_bar
        traj_flag(i)=2;
         t_norm_phi_bar=time(i);
        V(i)=uMax1*norm_phi+ c1* X(3:4,i)'*phi + c2*(uMax2*norm_v-V0_v)-V0_phi;  %*norm_v/uMax2*norm_v/uMax1
    elseif norm_phi>norm_phi_bar && norm_v<norm_v_bar
        traj_flag(i)=3;
        t_norm_phi_bar=time(i);
        V(i)=uMax1*norm_phi+ c1* X(3:4,i)'*phi + c2/(3-alpha2)*norm_v^(3-alpha2)-V0_phi; %*norm_v/uMax1
    else
        traj_flag(i)=4;
        V(i)=1/alpha_prime*norm(phi)^alpha_prime+ c1* X(3:4,i)'*phi + c2/(3-alpha2)*norm_v^(3-alpha2);
     end
     V(i)=1/alpha_prime*norm(phi)^alpha_prime+ c1* X(3:4,i)'*phi + c2/(3-alpha2)*norm_v^(3-alpha2);

    if norm_r>norm_r_bar2
       t_norm_r_bar2=time(i); 
    end
    
     if norm_r>norm_r_bar
       t_norm_r_bar=time(i); 
    end
    
    f=[X(3:4,i);U(:,i)-C_d*norm_v*X(3:4,i)];
    if i>1
    V_dot(i)=(V(i)-V(i-1))/dt;
    end
    V_dot_bar(i)=-k2*norm_v^(alpha2+1);
    time(i+1)=time(i)+dt;
%     if i>1
%     if V(i)-V(i-1)<1e-10 && V(i)<1e-5
%         break;
%     end
%     end
    if norm_r<1e-6 && norm_v<1e-6
        break;
    end
end
V(i+1)=V(i);
V_dot(i+1)=V_dot(i);
V_dot_bar(i+1)=V_dot_bar(i);
converge_time(j)=time(i);
converge_time_norm_r_bar(j)=t_norm_r_bar;
converge_time_norm_phi_bar(j)=t_norm_phi_bar;
%Upper bound on convergence time based on discontinuous control (first activate only velocity term until velocity to zero and the turn the position term on)
d1=norm_v0^(2-alpha2)/(k2*(2-alpha2));
t1=norm_v0^(1-alpha2)/(k2*(1-alpha2));
d2=sqrt(d1^2+norm_r0^2+2*X0(1:2)'*X0(3:4)*d1/norm_v0);
ftime=@(x) +norm_phi_bar-d2+(norm_v_bar_alpha/m1*(1-exp(m1*x))+norm_v_bar_alpha*x)-1/(k2*(2-alpha2))*abs(norm_v_bar_alpha*(1-exp(m1*x)))^(1-alpha2);
t2=fsolve(ftime,(d2-norm_r_bar-norm_v_bar_alpha/m1)/norm_v_bar_alpha); %(d2-norm_r_bar-norm_v_bar_alpha/m1)/norm_v_bar_alpha

ind0=find(rli<d2-norm_r_bar);
ind=max(ind0);

%     ftime=@(x) +norm_phi_bar-d2+rli(ind)+((Delta(ind)+vli(ind))/m_vl0(ind)*(exp(m_vl0(ind)*x)-1)-Delta(ind)*x)-1/(k2*(2-alpha2))*abs((Delta(ind)+vli(ind))*(exp(m_vl0(ind)*x))-Delta(ind))^(1-alpha2);
% t2=tli(ind)+fsolve(ftime,(d2-norm_r_bar-norm_v_bar_alpha/m1)/norm_v_bar_alpha); %(d2-norm_r_bar-norm_v_bar_alpha/m1)/norm_v_bar_alpha


%From finite time expression in unsaturated region
norm_v2=norm_v_bar_alpha*(1-exp(m1*t2));
norm_v2=abs(norm_v2);
V2=1/alpha_prime*norm_r_bar^(alpha_prime)+c2/(3-alpha2)*norm_v2^(3-alpha2);
t3=1/(c_bar*(1-2/(3-alpha2)))*V2^((1-2/(3-alpha2)));
converge_time_bar1(j)=t1+t2;%+t3;

%Upper bound on convergence time based on discontinuous control and velocity always pointing radially outward
m2=-k2*norm_v0^(alpha2-1);
t1=1/m2*log(-uMax1/(-uMax1+m2*norm_v0));
d1=uMax1/m2*(t1-exp(m2*t1)/m2)+norm_v0*exp(m2*t1)/m2-norm_v0/m2+uMax1/m2^2;
d2=abs(d1)+norm_r0;
ftime=@(x) +norm_phi_bar-d2+(norm_v_bar_alpha/m1*(1-exp(m1*x))+norm_v_bar_alpha*x)-1/(k2*(2-alpha2))*abs(norm_v_bar_alpha*(1-exp(m1*x)))^(1-alpha2);
t2=fsolve(ftime,(d2-norm_r_bar-norm_v_bar_alpha/m1)/norm_v_bar_alpha);

ind0=find(rli<d2-norm_r_bar);
ind=max(ind0);

%     ftime=@(x) +norm_phi_bar-d2+rli(ind)+((Delta(ind)+vli(ind))/m_vl0(ind)*(exp(m_vl0(ind)*x)-1)-Delta(ind)*x)-1/(k2*(2-alpha2))*abs((Delta(ind)+vli(ind))*(exp(m_vl0(ind)*x))-Delta(ind))^(1-alpha2);
% t2=tli(ind)+fsolve(ftime,(d2-norm_r_bar-norm_v_bar_alpha/m1)/norm_v_bar_alpha); %(d2-norm_r_bar-norm_v_bar_alpha/m1)/norm_v_bar_alpha



converge_time_bar2(j)=t1+t2;
end

converge_time_bound=1;
colors={[1,0,0],[0,0,0],[0,0,1]};
markers={'--','d','*','.'};
circCos=cos(0:pi/100:2*pi);
circSin=sin(0:pi/100:2*pi);

figure('units','normalized','outerposition',[.1 0.1 .8 .8]);
R0=norm_r_bar
plot(R0*circCos,R0*circSin)
hold on
R0=norm_phi_bar
plot(R0*circCos,R0*circSin)
hold on
R0=norm_r_bar2;
plot(R0*circCos,R0*circSin)
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
 plot(time(ind),-V(ind),markers{k},'color',colors{3},'markersize',3)

end
xlabel('t')
ylabel('V')
legend('$\dot{V}$','$\bar{\dot{V}}$','-V')

figure
plot(time,V_dot_bar)


figure
plot(Theta,converge_time,'b',Theta,converge_time_norm_phi_bar,'r',Theta,converge_time_bar1,'k',Theta,converge_time_bar2,'m')
xlabel('t')
ylabel('Convergence Time')


%%
%Fintie time controller for DI dynamics with quadratic drag term
%controller u = sat_u1(-k*r*||r||^(alpha1-1)) + sat_u2(-k*v*||v||^(alpha2-1)+C_d*||v||v)
%first slows down the agent to zero velocity and then starts the whole
%controller (but still the time of convergence remains to be determined)


uMax=10;
uMax1=7;
uMax2=uMax-uMax1;
k1=1;

alpha2=0.7;
alpha1=alpha2/(2-alpha2);
norm_r_bar = (uMax1/k1)^(1/alpha1);

V0=k1*alpha1/(1+alpha1)*(uMax1/k1)^((1+alpha1)/alpha1);
global C_d dt
C_d=0.2;

v_max=sqrt(uMax/C_d);

k2=(uMax-uMax2)/v_max^(alpha2);

f=@(x) C_d*x^2-k2*x^alpha2-uMax2;
norm_v_bar = abs(fsolve(f,sqrt(uMax2/C_d)));

dt=0.01;
NIter=50000;
NTheta=2;
for j=1:NTheta
    theta=(j-1)*2*pi/(NTheta-1);
    Theta(j)=theta;

X(:,1)=[148*rand(2,1);1;0];

flag_v=0;
time(1)=0;
for i=1:NIter
    norm_r=norm(X(1:2,i));
    if norm_r>1e-10
        U10=-k1*(X(1:2,i))*norm(X(1:2,i))^(alpha1-1);
    else
        U10=[0,0]';
    end
    norm_v=norm(X(3:4,i));
    norm_v_arr(i)=norm_v;
    if norm_v>1e-10 
        U20=-k2*(X(3:4,i))*norm(X(3:4,i))^(alpha2-1)+C_d*norm(X(3:4,i))*X(3:4,i);
    else
        U20=[0,0]';
    end
    norm_U10=norm(U10);
    norm_U20=norm(U20);
    if norm_v>1e-5 && flag_v==0
    t_v0=time(i);
    else
    flag_v=1;
    end
    uMax1_prime=uMax1*(min(1,time(i)-t_v0));
    if norm_U10<1e-10
        U(:,i) =  min(uMax2,norm_U20)*U20/norm_U20;
    elseif norm_U20<1e-10
        U(:,i) = min(uMax1_prime,norm_U10)*U10/norm_U10;
    else
    U(:,i) = min(uMax1_prime,norm_U10)*U10/norm_U10 + min(uMax2,norm_U20)*U20/norm_U20;
    end
    norm_U_arr(i)=norm(U(:,i));
    X(:,i+1)=modifiedDIDynamics(X(:,i),U(:,i));
    
    %Find the Lyapunov function at time i
    
    if norm_r<norm_r_bar
        V(i)=k1/(1+alpha1)*norm_r^(1+alpha1) + 0.5*norm_v^2;
        dV_by_dX=[k1*X(1:2,i)*norm_r^(alpha1-1);X(3:4,i)];
    else
        V(i)=uMax1*norm_r + 0.5*norm_v^2 -V0;
        dV_by_dX=[uMax*X(1:2,i)*norm_r;X(3:4,i)];
    end
    
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
    
    f=[X(3:4,i);U(:,i)-C_d*norm_v*X(3:4,i)];
    V_dot(i)=dV_by_dX'*f;
    V_dot_bar(i)=-k2*norm_v^(alpha2+1);
    time(i+1)=time(i)+dt;
%     if i>1
%     if V(i)-V(i-1)<1e-10 && V(i)<1e-5
%         break;
%     end
%     end
    if norm_r<1e-6 && norm_v<1e-6
        break;
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

figure
plot(time,X(1,:),time,X(2,:))

figure
%subplot(2,1,1)
plot(time,V)
xlabel('t')
ylabel('V')

figure
plot(time,V_dot,time,-V)%,time,V1_dot,time,V2_dot)
xlabel('t')
ylabel('V')
legend('V_dot','-V')

figure
plot(time,V_dot_bar)


figure
plot(Theta,convergence_time)
xlabel('t')
ylabel('Convergence Time')

