%plot appollonius circles corresponding to the two defenders and one
% ------------------------------
% Author: Vishnu S. Chipade
% Email: vishnuc@umich.edu
% ------------------------------

%attacker

rA=[rD(1,2)/2,10]';
vA=1.9;
ang=[0,60]'*pi/180;
rho_AD=10;
rD=rA+rho_AD*[cos(ang(1)),sin(ang(1));cos(ang(2)),sin(ang(2));]';
vD=1;

w=vD/vA;

circCos=cos(0:pi/100:2*pi);
circSin=sin(0:pi/100:2*pi);

colors={[0,0,1],[0.5,0.1,.1]};

%Plot the points
figure
hold on
%Attacker
plot(rA(1),rA(2),'ro')

%find the appollonius circle and plot
for i=1:2
    %Plot the defenders
    plot(rD(1,i),rD(2,i),'o','color',colors{i})
    text(rD(1,i),rD(2,i),['$\mathcal{D}_',num2str(i),'$'])
    %find the circle parameters
    rC_App(:,i)=1/(1-w^2)*[rD(1,i)-w^2*rA(1),rD(2,i)-w^2*rA(2)]';
    R_App(i)=norm(rD(:,i)-rA)*w/(1-w^2);
    %plot circle
    plot(rD(1,i)+R_App(i)*circCos,rD(2,i)+R_App(i)*circSin,'color',colors{i});
end

