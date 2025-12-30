% ------------------------------
% Author: Vishnu S. Chipade
% Email: vishnuc@umich.edu
% ------------------------------

%SetPlotDefaults;
set(groot,'defaulttextinterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');


global options

options = optimset('Display','off','MaxIter',1000,'TolFun',1e-3);

global dt Rs  alpha rho_B rho_c_A rho_P rho_D rho_A;

dt=0.01;  %time step

%% General
Rs=1000;   %Sensing radius
rho_P=45;  %Radius of the protected area

%Initial Positions of the protected and safe area
global rP
rP=[0,0]';

global rS rho_S
rS=[300,500]';

%% Attackers

%Parameters for the attackers
rho_A=0.5;
rho_c_A=150;

%Control Limits
global v_maxA u_maxA C_d
C_d=0.2;
vma=6;
if C_d~=0
    uma=C_d*vma.^2;
    %temporary safe distances (not based on the quadratic drag term
    %analysis)
    rhoA_safe=log((uma+C_d*(2*vma)^2)/uma)/(2*C_d);     %additional aafety distance required between the attackers due to bounded acceleration
else
    uma=0.35*vma;
    rhoA_safe=(vma(1)+vma(1))^2/(2*uma(1));
end

%Parameters for the vector fields for the attacker corresponding to the
%other attackers
global R_m_AA R_bar_AA R_u_AA
global A_A_A B_A_A C_A_A D_A_A
R_bar_AA1=70*rho_A;
R_bar_AA2=80*rho_A;
R_m_AA=1*(rho_A+rho_A)+rhoA_safe+1;
R_bar_AA=R_m_AA+R_bar_AA1+rhoA_safe;
R_u_AA=R_m_AA+R_bar_AA2+rhoA_safe;
dR_AA_cube=(R_u_AA-R_bar_AA)^3;
A_A_A=2/dR_AA_cube;
B_A_A=-3*(R_u_AA+R_bar_AA)/dR_AA_cube;
C_A_A=6*R_u_AA*R_bar_AA/dR_AA_cube;
D_A_A=R_u_AA^2*(R_u_AA-3*R_bar_AA)/dR_AA_cube;


%Initialize the attackers

NA=3;  %number of attackers
v_maxA=vma*ones(1,NA);  %attacker velocity bound array
u_maxA=uma*ones(1,NA);  %attacker acceleration bound array


axS=[1,1]';
rDmin=rho_P+10;
rDmin=7;

global RA0 rho_Acon
RA0=R_m_AA+4;
global Rii00 Rik00 Rjk00
Rii00=RA0*sqrt(2*(1-cos(2*pi/NA)))*ones(1,NA);  %For formation potential
Rik00=5*Rii00;
dRA0=1;

%Initialize the attackers
rA0=[-750,280]';
%rA0=[-955,700]';
%rA0=[-955,140]';
%rA0=[-875,200]';
thetaA0=atan2(rA0(2)-rP(2),rA0(1)-rP(1))+pi/3;
rA(:,1)=rA0;
RA01=RA0;
% rA(:,2)=rA0+RA01*[cos(thetaA0-pi/10), sin(thetaA0-pi/10)]';
% rA(:,3)=rA0+RA01*[cos(thetaA0+2*pi/10), sin(thetaA0+2*pi/10)]';
%XA0(:,1:3)=[rA(:,1:3);zeros(2,3)];
for i=1:NA
    thetaA=thetaA0+2*pi/NA*(i-1);
    
    %RA=sqrt(2)*(RA0)+RA0*(i-4)+1*randn;
    RA=RA0;%*(i-1);%+1*randn;
    %RA=RA0;
    %RA=5*RA0*rand(1);
    rA(:,i)=rA0+[RA*cos(thetaA), RA*sin(thetaA);]';
    vA(:,i)=-(rA(:,i)-rP)*0.01*v_maxA(i)/norm(rA(:,i)-rP);%+rand(2,1);
    vA(:,i)=[0,0]';
    %vA(:,i)=v_maxA(i)*([-1,-1]'+2*rand(2,1));
    XA0(:,i)=[rA(:,i);vA(:,i)];
end
rA_follow=rA-rA(:,1);

% fread=fopen('/home/msr/msr_ws/src/rotors_simulator/rotors_gazebo/launch/test3.launch','r');
% fwrite=fopen('test.launch','w');
% 
% while ~feof(fread)  
%     fline=fgets(fread);
%     num=textscan(fline,'<group ns="attacker%d">');
%     while ~isempty(num{1}) && ~feof(fread) 
%         fprintf(fwrite,fline);
%          fline=fgets(fread);
%         num_temp=textscan(fline,'<group ns="attacker%d">');
%         while isempty(num_temp{1}) && ~feof(fread)           
%         num_x=textscan(fline,'<arg name="x" value="%f"/>');
%         num_y=textscan(fline,'<arg name="y" value="%f"/>');
%         num_waypoint=textscan(fline,'<node name="waypoint_publisher" pkg="rotors_gazebo" type="waypoint_publisher" output="screen" args="%f %f %f %f %f"/>');
%             if ~isempty(num_x{1})
%                 fprintf(fwrite,'      <arg name="x" value="%f"/>\n',rA(1,num{1}));
%             elseif ~isempty(num_y{1})
%                 fprintf(fwrite,'      <arg name="y" value="%f"/>\n',rA(2,num{1}));
%             elseif ~isempty(num_waypoint{1})
%                 fprintf(fwrite,'<node name="waypoint_publisher" pkg="rotors_gazebo" type="waypoint_publisher" output="screen" args="%f %f %f %f %f"/>\n',[rA(1,num{1}),rA(2,num{1}), num_waypoint{3:5}]);
%             else                
%                 fprintf(fwrite,fline);
%             end
%             fline=fgets(fread);
%         num_temp=textscan(fline,'<group ns="attacker%d">');
%         end 
%         num=num_temp;        
%     end
% end
% 
% 
% fclose(fread);
% 
% fclose(fwrite);

%Create entries for lanuch file


%XA0(:,NA)=[-700,0,0,0]';
rAcm0=sum(XA0(1:2,:),2)/NA;
vAcm0=sum(XA0(3:4,:),2)/NA;

thetaAcm0=atan2(rAcm0(2)-rP(2),rAcm0(1)-rP(1));
%parameters for the defenders formation
global rho_sn
rho_Acon=RA0+rhoA_safe+dRA0;  %radius of the connectivity region

rho_S=4*rho_Acon;   %radius of the safe area


%% Defenders
rho_D=0.5;
rho_sn_max=28;
ND=ceil(pi/acos(rho_Acon/rho_sn_max));  %number of defenders
N=NA+ND;

global v_maxD v_maxDC u_maxD u_maxD1 u_maxD2 u_maxDr1 u_maxDr2 rho_safe obs dthetai alphaD_v
global rhoAD_safe rhoD_safe
alphaD_v=0.5;
kDv=.05;
alphaD_r=alphaD_v/(2-alphaD_v);
umdr2=uma+kDv*vma;
umdr1=0.35*umdr2;
eta_max=vma+.2*vma;
umd2=max(umdr1+umdr2,C_d*eta_max^2+kDv*eta_max^(alphaD_v));
umd1=1.001*umd2;
umd=umd1+umd2;
if C_d~=0
    vmd=sqrt(umd/C_d);
    rhoD_safe=log((umd+C_d*(2*vmd)^2)/umd)/(2*C_d);   %additional aafety distance required between the defenders due to bounded acceleration
    rhoAD_safe=log((uma+C_d*(vma+vmd)^2)/uma)/(2*C_d);  %additional safty for the attackers to stay safe from the defenders due to bounded acceleration
    rho_safe = rhoD_safe;
    vmdc=sqrt((sqrt(umd^2/(1/rho_safe^2+C_d^2))));
else
    u_maxD=0.4*vmd;
    rhoD_safe=(vmd+vmd)^2/(2*umd);
    rhoAD_safe=(vma+vmd)^2/(2*uma);
    rho_safe = rhoD_safe;
    vmdc=sqrt(umd*rhoD_safe);
end
v_maxD=vmd*ones(1,ND);
v_maxDC=vmdc*ones(1,ND);
u_maxD=umd*ones(1,ND);
u_maxD1=umd1*ones(1,ND);
u_maxD2=umd2*ones(1,ND);
u_maxDr1=umdr1*ones(1,ND);
u_maxDr2=umdr2*ones(1,ND);
obs.rho_safe=rho_safe;
dthetai=acos(1-(2*rho_D)^2/(2*rho_safe^2));   %angular shift for two agents colliding on a ciruclar arc segment

Rjk00=Rik00(1)*ones(1,ND);


global largeP
largeP=100000;
%Parameters for finite time control
global kr0 kr1 kr2 kv1 kv2 RD_con
%kr0=v_maxD(1)-v_maxA(1)-0.5;
kr2=.5;
f0=@(x) (1-tanh(x)^2)-kr2*tanh(x)/x;
RD_con=fsolve(f0,1.5,options);
%kr1=kr0*tanh(RD_con)/RD_con^kr2;
kr1=0.5;
kv1=0.1;
kv2=0.5;

%control gains for attackers
global kADr kAFr kAOr kAOr2 kAPr kAPv kAFv alphaAFv kAOv kAOv2 alphaAOv kADv alphaADv
kAFr=10;
kAOr=1.5;
kAOr2=0.5;
kADr=.5; %kADr=1;
kADv=.5;  %kADv=2;
kAFv=.1;
alphaAFv=1;
kAOv=2;
kAOv2=.2;
alphaAOv=0.5;
alphaADv=0.5;
kAPr=0.003;
kAPv=0.1;

%finite time gains for the defenders
global kDOr kDDr kDOv1 kDOv2 alphaDOv kDDv alphaDDv
global kDFr  kDFr2 kDFv alphaDFr alphaDFv kDRr kDRv
global kDFphi kDFphiv kDFphir kDFphid
global kDDesr kDDesv
kDOr=1;
kDDr=.5;
kDOv1=2;
kDOv2=0.08;
kDFr=5.5;
kDFr2=.2;
kDFv=.45;
alphaDFv=.9;
alphaDFr=alphaDFv/(2-alphaDFv);
alphaDOv=1;
kDDv=0.3;
alphaDDv=0.5;
kDRr=.0001;
kDRv=0.1;
kDFphi=0.003;
kDFphid=0.007;
kDFphir=0.02;
kDFphiv=0.1;
kDDesr=1;
kDDesv=1;


%bound on the convergence error during tracking
A_tilde=[zeros(2),eye(2);-kDDesr*eye(2),-kDDesv*eye(2)];
Qd=eye(4);
Pd=lyap(A_tilde,Qd);
c1=min(eig(Pd));
c2=max(eig(Pd));
c3=min(eig(Qd));
c4=2*max(eig(Pd));
p=c4/c3*sqrt(c2/c1);
bd=15;

%Convex polygonal Obstalces
% global rO
% rO=[];
% rO=10*[9,22,2,3;-6,18,3,4;11,5,2,2;15,43,3,3;-2,45,3,4;12,60,4,3]'; %center (xO,yO), width (wO) and height (hO) of the rectangular obstacles
% rO=[110,350,100,80;-245,455,110,210;95,-210,50,140;-165,-110, 130, 100]';  %;-160,258,50,41
%rVO={[0,1;6,1;6,5;0,5;]',[14,3;20,3;20,9;14,9]',[4,-14;11,-14;11,-8;4,-8;]',[7,16;14,16;14,22;7,22;]'};
rVO={[60,310;160,310;160,390;60,390;]',[-300,350;-190,350;-190,560;-300,560]',[70,-280;120,-280;120,-140;70,-140;]',[-230,-160;-100,-160;-100,-60;-230,-60;]',...
    [-600,400;-500,400;-500,550;-600,550]', [-700,0;-570,0;-570,90;-700,90;]',[-400,-500;-250,-500;-250,-400;-400,-400;]',[-500,850;-300,850;-300,970;-500,970;]',...
    [400,100;650,100;650,260;600,260;]'};
rVO={[100,410;200,410;200,550;100,550;]',[-350,450;-240,450;-240,660;-350,660]',[-230,-160;-100,-160;-100,-60;-230,-60;]',...
    [-800,400;-600,400;-700,550;-800,550]', [-650,0;-520,0;-520,90;-620,150 ;-680,90;]',[300,-50;550,-50;500,160;350,160;]'};


%Rectangular obstacles
rCO=[150, 388; -310, 388; -179, -120; -650, 300; -450, 77; 340, 77.5  ]';
lxO=[65,65,65,130,65,45];
lyO=[85,85,45,45,85,65];

for k=1:6
    rVO{k}=[rCO(1,k)-lxO(k)/2,rCO(2,k)-lyO(k)/2;  rCO(1,k)+lxO(k)/2,rCO(2,k)-lyO(k)/2;   rCO(1,k)+lxO(k)/2,rCO(2,k)+lyO(k)/2;   rCO(1,k)-lxO(k)/2,rCO(2,k)+lyO(k)/2]';
end
%[-400,-500;-250,-500;-250,-400;-400,-400;]',[-500,850;-300,850;-300,970;-500,970;]'

obs.rVO=rVO;

if (1)
    figure
    fontSize=18;
    for k=1:size(rVO,2)
        rVO1{k}=[rVO{k},rVO{k}(:,1)];
        XO{k}= rVO1{k}(1,:);
        YO{k}= rVO1{k}(2,:);
        polyin=polyshape(rVO{k}(1,:),rVO{k}(2,:));
        [rCO2(1,k),rCO2(2,k)]=centroid(polyin);
        hold on;
        fill(XO{k},YO{k},[0.2 0.2 0.2])
        text(sum(rVO{k}(1,:))/length(rVO{k}(1,:)),sum(rVO{k}(2,:))/length(rVO{k}(2,:)),['$\mathcal{O}_',num2str(k),'$'],'fontsize',0.8*fontSize,'color',[1,1,1]);
    end
end
obs.rVO1=rVO1;
NO=length(rVO);
obs.NO=NO;
obs.rCO2=rCO2;

%%
%For formation orientation
global R_bar_AcOc R_u_AcOc A_Ac_Oc B_Ac_Oc C_Ac_Oc D_Ac_Oc
for k=1:NO
    NVOk=length(rVO{k}(1,:));
    rOc=sum(rVO{k},2)/NVOk;
    for j=1:NVOk
        dVO(j)=norm(rOc-rVO{k}(:,j));
    end
    R_bar_AcOc(k)=rho_Acon+max(dVO)+5;
    R_u_AcOc(k)=rho_Acon+max(dVO)+155;
    dR_AcOc_cube=(R_u_AcOc(k)-R_bar_AcOc(k))^3;
    A_Ac_Oc(k)=2/dR_AcOc_cube;
    B_Ac_Oc(k)=-3*(R_u_AcOc(k)+R_bar_AcOc(k))/dR_AcOc_cube;
    C_Ac_Oc(k)=6*R_u_AcOc(k)*R_bar_AcOc(k)/dR_AcOc_cube;
    D_Ac_Oc(k)=R_u_AcOc(k)^2*(R_u_AcOc(k)-3*R_bar_AcOc(k))/dR_AcOc_cube;
end

%Parameters for the vector fields for the attacker corresponding to the
%defenders
global R_m_AD R_bar_AD R_u_AD
global A_A_D B_A_D C_A_D D_A_D
R_bar_AD1=40*rho_A;  %65
R_bar_AD2=50*rho_A;  %75
R_m_AD=1*(rho_A+rho_D)+rhoAD_safe+2;
%R_m_AD=0;  % For testing purposes
R_bar_AD=R_m_AD+R_bar_AD1+rhoAD_safe;
R_u_AD=R_m_AD+R_bar_AD2+rhoAD_safe;
dR_AD_cube=(R_u_AD-R_bar_AD)^3;
A_A_D=2/dR_AD_cube;
B_A_D=-3*(R_u_AD+R_bar_AD)/dR_AD_cube;
C_A_D=6*R_u_AD*R_bar_AD/dR_AD_cube;
D_A_D=R_u_AD^2*(R_u_AD-3*R_bar_AD)/dR_AD_cube;

global R_m_AD2 R_bar_AD2 R_u_AD2
global A_A_D2 B_A_D2 C_A_D2 D_A_D2
R_m_AD2=11*R_m_AD;%1*(rho_A+rho_D)+rhoAD_safe;
%R_m_AD=0;  % For testing purposes
R_bar_AD2=R_m_AD2+R_bar_AD1+rhoAD_safe;
R_u_AD2=R_m_AD2+R_bar_AD2+rhoAD_safe;
dR_AD_cube2=(R_u_AD2-R_bar_AD2)^3;
A_A_D2=2/dR_AD_cube2;
B_A_D2=-3*(R_u_AD2+R_bar_AD2)/dR_AD_cube2;
C_A_D2=6*R_u_AD2*R_bar_AD2/dR_AD_cube2;
D_A_D2=R_u_AD2^2*(R_u_AD2-3*R_bar_AD2)/dR_AD_cube2;

global Rij0 Rjj0
Rij0=R_u_AD*ones(1,ND);
Rjj0=20*ones(1,ND);



%%

%RAD_max=(1-rad0)*R_m_AD+rad0*R_bar_AD; %The locations of the defenders placed around the attacker
RAD_max=0.2*R_bar_AD+.8*R_m_AD;
rho_Fmax=RAD_max+rho_D;

%Parameters for the vector fields guiding the entire formation
alpha=1;  %Coefficient multiplied in finding n for superelliptic contours
global nO A B C D aO bO E_bar_O E_m_O E_u_O GO
GO=500;
R_safe=5*rho_D;
E_OF1=2;  %repulsive
E_OF2=16;  %blending of repulsive and attractive

for k=1:NO
    w(k)=rVO{k}(1,2)-rVO{k}(1,1);
    h(k)=rVO{k}(2,3)-rVO{k}(2,2);
    w_bar(k)=w(k)+2*(rho_Fmax+R_safe);
    h_bar(k)=h(k)+2*(rho_Fmax+R_safe);
    f1=@(x) x-0.5*(w_bar(k)/w(k))^(2/(1-exp(-x)))-0.5*(h_bar(k)/h(k))^(2/(1-exp(-x)))+1;
    E_m_O(k)=fsolve(f1,rho_Fmax+sqrt(w(k)^2+h(k)^2)/2,options);
    nO(k)=1/(1-exp(-alpha*E_m_O(k)));   %power of super-elliptic contour
    E_u_O(k)=E_m_O(k)+E_OF2; %Maximum elliptic distance of effective circulant vector field
    
    E_bar_O(k)=E_m_O(k)+E_OF1;  %Minimum elliptic distance after which two vector fields are combined
    dREOcube=(E_u_O(k)-E_bar_O(k))^3;
    A(k)=2/dREOcube;
    B(k)=-3*(E_u_O(k)+E_bar_O(k))/dREOcube;
    C(k)=6*E_u_O(k)*E_bar_O(k)/dREOcube;
    D(k)=E_u_O(k)^2*(E_u_O(k)-3*E_bar_O(k))/dREOcube;
    aO(k)=w(k)/2*(2^(1/(2*nO(k))));
    bO(k)=h(k)/2*(2^(1/(2*nO(k))));
end

%Parameters for the vector fields for the attacker corresponding to the
%obstacles
global R_m_AO R_bar_AO R_u_AO R_v_AO
global A_A_O B_A_O C_A_O D_A_O
global A_bar_A_O B_bar_A_O C_bar_A_O D_bar_A_O
R_bar_O1=7*rho_A;  %repulsive
R_bar_O2=17*rho_A;  %blending of repulsive and attractive
R_bar_O3=12;
for j=1:NA
    for k=1:NO
        %f1=@(x) x-0.5*((w(k)+2*rho_A))^(2/(1-exp(-x)))-0.5*((h(k)+2*rho_A)/h(k))^(2/(1-exp(-x)))+1;
        rho_bar=2*10*rho_A;
        R_m_AO(j,k)=0;
        R_bar_AO(j,k)=2*30*rho_A;
        R_u_AO(j,k)=2*40*rho_A;
        R_v_AO(j,k)=2*40*rho_A;
        dR_AO_cube=(R_u_AO(j,k)-R_bar_AO(j,k))^3;
        A_A_O(j,k)=2/dR_AO_cube;
        B_A_O(j,k)=-3*(R_u_AO(j,k)+R_bar_AO(j,k))/dR_AO_cube;
        C_A_O(j,k)=6*R_u_AO(j,k)*R_bar_AO(j,k)/dR_AO_cube;
        D_A_O(j,k)=R_u_AO(j,k)^2*(R_u_AO(j,k)-3*R_bar_AO(j,k))/dR_AO_cube;
        
        dR_bar_AO_cube=(R_v_AO(j,k)-R_u_AO(j,k))^3;
        A_bar_A_O(j,k)=2/dR_bar_AO_cube;
        B_bar_A_O(j,k)=-3*(R_v_AO(j,k)+R_u_AO(j,k))/dR_bar_AO_cube;
        C_bar_A_O(j,k)=6*R_v_AO(j,k)*R_u_AO(j,k)/dR_bar_AO_cube;
        D_bar_A_O(j,k)=R_v_AO(j,k)^2*(R_v_AO(j,k)-3*R_u_AO(j,k))/dR_bar_AO_cube;
    end
end



%Parameters for the vector fields for the attacker corresponding to the
%obstacles
if(0)
    global R_m_AO R_bar_AO R_u_AO
    global A_A_O B_A_O C_A_O D_A_O
    R_bar_O1=15*rho_A;  %repulsive
    R_bar_O2=20*rho_A;  %blending of repulsive and attractive
    for k=1:NO
        R_m_AO(k)=1/2*sqrt(w_bar(k)^2+h_bar(k)^2);
        R_bar_AO(k)=R_m_AO(k)+R_bar_O1;
        R_u_AO(k)=R_bar_AO(k)+R_bar_O2;
        dR_AO_cube=(R_u_AO(k)-R_bar_AO(k))^3;
        A_A_O(k)=2/dR_AO_cube;
        B_A_O(k)=-3*(R_u_AO(k)+R_bar_AO(k))/dR_AO_cube;
        C_A_O(k)=6*R_u_AO(k)*R_bar_AO(k)/dR_AO_cube;
        D_A_O(k)=R_u_AO(k)^2*(R_u_AO(k)-3*R_bar_AO(k))/dR_AO_cube;
    end
end
%Parameters for the vector fields for the defenders
global R_m_DD R_m_DDO R_bar_DD R_u_DD
global A_D_O B_D_O C_D_O D_D_O A_D_D B_D_D C_D_D D_D_D
global R_m_DO R_bar_DO R_u_DO

E_bar_O1=10*rho_D;  %repulsive
E_bar_O2=20*rho_D;  %blending of repulsive and attractive
E_bar_O3=5;
for j=1:ND
    for k=1:NO
        rho_bar=2*10*rho_D;
        R_m_DO(j,k)=0;
        R_bar_DO(j,k)=2*30*rho_D;
        R_u_DO(j,k)=2*40*rho_D;
        dR_DO_cube=(R_u_DO(j,k)-R_bar_DO(j,k))^3;
        A_D_O(j,k)=2/dR_DO_cube;
        B_D_O(j,k)=-3*(R_u_DO(j,k)+R_bar_DO(j,k))/dR_DO_cube;
        C_D_O(j,k)=6*R_u_DO(j,k)*R_bar_DO(j,k)/dR_DO_cube;
        D_D_O(j,k)=R_u_DO(j,k)^2*(R_u_DO(j,k)-3*R_bar_DO(j,k))/dR_DO_cube;
    end
end
%For the formation of defenders
if (0)
    global E_m_DcO E_bar_DcO E_u_DcO E_v_DcO
    global A_Dc_O B_Dc_O C_Dc_O D_Dc_O
    global A_bar_Dc_O B_bar_Dc_O C_bar_Dc_O D_bar_Dc_O
    for k=1:NO
        rho_bar=2*rho_sn;
        E_m_DcO(1,k)=0.5*((w(k)+rho_bar)/w(k))^(2*nO(k))+0.5*((h(k)+rho_bar)/h(k))^(2*nO(k))-1;
        rho_bar=2*1.05*rho_sn;
        E_bar_DcO(1,k)=0.5*((w(k)+rho_bar)/w(k))^(2*nO(k))+0.5*((h(k)+rho_bar)/h(k))^(2*nO(k))-1;
        rho_bar=2*1.1*rho_sn;
        E_u_DcO(1,k)=0.5*((w(k)+rho_bar)/w(k))^(2*nO(k))+0.5*((h(k)+rho_bar)/h(k))^(2*nO(k))-1;
        rho_bar=2*1.4*rho_sn;
        E_v_DcO(1,k)=0.5*((w(k)+rho_bar)/w(k))^(2*nO(k))+0.5*((h(k)+rho_bar)/h(k))^(2*nO(k))-1;
        dE_DcO_cube=(E_u_DcO(1,k)-E_bar_DcO(1,k))^3;
        A_Dc_O(1,k)=2/dE_DcO_cube;
        B_Dc_O(1,k)=-3*(E_u_DcO(1,k)+E_bar_DcO(1,k))/dE_DcO_cube;
        C_Dc_O(1,k)=6*E_u_DcO(1,k)*E_bar_DcO(1,k)/dE_DcO_cube;
        D_Dc_O(1,k)=E_u_DcO(1,k)^2*(E_u_DcO(1,k)-3*E_bar_DcO(1,k))/dE_DcO_cube;
        
        dE_bar_DcO_cube=(E_v_DcO(1,k)-E_u_DcO(1,k))^3;
        A_bar_Dc_O(1,k)=2/dE_bar_DcO_cube;
        B_bar_Dc_O(1,k)=-3*(E_v_DcO(1,k)+E_u_DcO(1,k))/dE_bar_DcO_cube;
        C_bar_Dc_O(1,k)=6*E_v_DcO(1,k)*E_u_DcO(1,k)/dE_bar_DcO_cube;
        D_bar_Dc_O(1,k)=E_v_DcO(1,k)^2*(E_v_DcO(1,k)-3*E_u_DcO(1,k))/dE_bar_DcO_cube;
    end
end
%%
%for other defenders
global R_m2_DD
R_bar_DD1=6*rho_D;
R_bar_DD2=6*rho_D;
R_m_DD=1*(rho_D+rho_D)+rhoD_safe+2;
R_m_DDO=3*R_m_DD;
R_m2_DD=100;
R_bar_DD=R_m_DD+R_bar_DD1+rhoD_safe;
R_u_DD=R_bar_DD+R_bar_DD2+rhoD_safe;
dR_DD_cube=(R_u_DD-R_bar_DD)^3;
A_D_D=2/dR_DD_cube;
B_D_D=-3*(R_u_DD+R_bar_DD)/dR_DD_cube;
C_D_D=6*R_u_DD*R_bar_DD/dR_DD_cube;
D_D_D=R_u_DD^2*(R_u_DD-3*R_bar_DD)/dR_DD_cube;

rho_sn=rho_Acon+R_m_DD+6;   %radius of the stringNet




%%

%Initialize the defenders
global rSD_goal
%rD0=[-120,-220]';
rD0=[-150,-300]';
for j=1:ND
    RD0=20+(j-1)*8;%3*rDmin;
    thetaD=(j-1)*(2*pi-pi/3)/ND;
    rD(:,j)=rD0+RD0*[cos(thetaD), sin(thetaD);]';
    XD0(:,j)=[rD(:,j);0;0];
    rSD_goal(:,j)=rS+RD0*[cos(thetaD), sin(thetaD);]';
end

rD=[-195,-195;-60,-165;-500,0;300,150;]';%-150,500]';
%rD=[-190,-250;-180,-258;-800,600;270,450; -150,500]';
%for direct seeking phase
%rD=[-400.8,253;-352,280;-454,268;-337,333;-364,381]';
XD0=[rD;zeros(2,ND)];

%add virtual defender at rDFcm

%%
%Define the system parameters on pelican quadrotors;
m=1; %mass
g=9.81;
%Inertia
I_xx=0.01;
I_yy=0.01;
I_zz=0.02;
I=diag([I_xx,I_yy,I_zz]);
inv_I=inv(I);
a_l=0.21;  %arm length of the quadrotor
K_f=9.9865*10^(-6);
K_m=1.6*10^-2;
alloc_mat= [K_f,  K_f,  K_f,  K_f;...
    0,  a_l*K_f,  0,  -a_l*K_f;...
    -a_l*K_f,  0,  a_l*K_f,  0;...
    K_f*K_m,  -K_f*K_m, K_f*K_m, -K_f*K_m];  % allocation matrix
inv_alloc_mat=inv(alloc_mat);
acc2omega=inv_alloc_mat*[1,zeros(1,3); zeros(3,1), I];

%Low level controller Gain Parameters
pos_gain=.9*diag([4,4,5]);
vel_gain=diag([2.7,2.7,2.7]);
att_gain=diag([.6,.6,0.035]);
%att_gain=diag([1,1,0.035]);
ang_rate_gain=diag([0.22,0.22,0.05]);
att_gain_mar=inv_I*diag([1,1,0.035]);
ang_rate_gain_mar=inv_I*diag([0.22,0.22,0.05]);

global params;
params.g=9.81;
params.m=m;
params.I=I;
params.omega_max=834;
params.alloc_mat=alloc_mat;
params.inv_alloc_mat=inv_alloc_mat;
params.acc2omega=acc2omega;
params.pos_gain=pos_gain;
params.vel_gain=vel_gain;
params.att_gain=att_gain;
params.ang_rate_gain=ang_rate_gain;
params.att_gain_mar=att_gain_mar;
params.ang_rate_gain_mar=ang_rate_gain_mar;
