%List of testings
% ------------------------------
% Author: Vishnu S. Chipade
% Email: vishnuc@umich.edu
% ------------------------------

%1) variance plotting
%2) Potential function and it's gradient

%Test Variance
N=10;
r=rand(N,2)
rcm=sum(r,1)/N;
sigma2=0;
sigma=0;
maxR=-Inf;
rdist=pdist(r);
for i=1:N
    d=norm((r(i,:)-rcm));
    sigma2=sigma2+d^2;
    sigma=sigma+d;
    dist(i)=d;
    if d>maxR
        maxR=d;
        indMax=i;
    end
end
medianR=median(dist);
sigma2=sigma2/N;
sigma=sigma/N;
figure
plot(r(:,1),r(:,2),'o',rcm(1),rcm(2),'m*', [rcm(1),r(indMax,1)],[rcm(2),r(indMax,2)],'--')
hold on;
Rad=sqrt(sigma2);
plot(rcm(1)+Rad*cos(0:pi/50:2*pi),rcm(2)+Rad*sin(0:pi/50:2*pi),'r')
Rad=sigma;
plot(rcm(1)+Rad*cos(0:pi/50:2*pi),rcm(2)+Rad*sin(0:pi/50:2*pi),'r--')
Rad=medianR;
plot(rcm(1)+Rad*cos(0:pi/50:2*pi),rcm(2)+Rad*sin(0:pi/50:2*pi),'m--')
%%

%Test Potential function and gradient
x=0:0.1:100;
alpha=0.5;
c=0.1;
r0=10;
rm=1;
V=log((.5*x-r0)./(x-rm)+(x-rm)./(.5*x-r0));
V1=V.^(alpha);
V2=V.^(1+alpha);
norm_nabla_V=abs(((x-rm).^2-r0^2)./((x-rm).^2+r0^2));
norm_nabla_V1=norm_nabla_V.^(alpha-1);
V2=log(r0./abs(x-rm)+abs(x-rm)./r0);
nabla_V2=((x-rm).^2-r0^2)./((x-rm).^2+r0^2)./abs(x-rm)
figure;
plot(x,V,'b-',x,norm_nabla_V,'b--')

if (0)
    figure
plot(x,V2,x,nabla_V2);
    figure;
    plot(x,V1,'b',x,norm_nabla_V1,'b--');
    hold on;
    plot(x,V1-norm_nabla_V1,'r',x,-c*V1,'r--',x,-c*V2,'m--')
end

%with absolute input argument
V=log(r0./x+x./r0);

%%
%%Saturated input

X(:,1)=[15,5]';
uMax=10000;
dt=0.01;
time(1)=0;
NIter=5000;
for i=1:NIter
    g=min(uMax,abs(-X(2,i)-X(1,i)))/abs(-X(2,i)-X(1,i));
    f=[X(2,i),g*(-X(2,i)-X(1,i))]';
    X(:,i+1)=X(:,i)+dt*f;
    V(i)=0.5*sum(X(:,i).^2);
    V_dot(i)=X(:,i)'*f;
    V1_dot(i)=(1-g)*X(2,i)*X(1,i);
    V2_dot(i)=-g*X(2,i)^2;
    time(i+1)=time(i)+dt;
end
V(i+1)=V(i);
V_dot(i+1)=V_dot(i);
V1_dot(i+1)=V1_dot(i);
V2_dot(i+1)=V2_dot(i);
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


%% Bound on error in double integrator system with asymptotic/exponential convergence
k1=5; k2=5;
A=[0,1;-k1,-k2];
theta=0.99999;
delta=0.3;
Q=2*eye(2);
eigQ=eig(Q)
P=lyap(A,Q)
eigP=eig(P)
C1=min(eigP);
C2=max(eigP);
C3=min(eigQ);
C4=2*C2;
bf=C4/C3*sqrt(C2/C1);
bf1=1/bf;
r=delta/bf1;
b=bf*delta/theta


%%
% Plot logarithmic potential and its gradient


%%
%Shortest path in presence of obstacle
rO=[5,0]';
RO=3;
rP1=[0,0]';
rP2=[10,0]';



%%
%Move along a superelliptic curve
k=1;
rO=[0;0];
w=7;
h=5;
n=1;
rho_A=1;
a=w/2*(2^(1/(2*n)));
b=h/2*(2^(1/(2*n)));
RcO=0.5*sqrt(w^2+h^2);
Acirc=pi*RcO^2;
v0=10;
EO=(w/2/a)^(2*n)+(h/2/b)^(2*n)-1+5;
x0=a*(EO+1)^(1/(2*n));
ef=@(x) (EO+1-(x./a).^(2*n)).^(1/2/n);
Aell=4*b*integral(ef,0,x0);
X(:,1)=[x0,0,0,v0]';
xsup=a*(EO+1)^(1/(2*n))*cos(0:pi/50:pi/2).^(1/n);
ysup=b*(EO+1)^(1/(2*n))*sin(0:pi/50:pi/2).^(1/n);
xsup=[xsup,-flip(xsup),-xsup,flip(xsup)];
ysup=[ysup,flip(ysup),-ysup,-flip(ysup)];
xsup=xsup+rO(1,k);
ysup=ysup+rO(2,k);
circCos=cos(0:pi/100:2*pi);
circSin=sin(0:pi/100:2*pi);

NIter=5000;
dt=0.001;
time(1)=0;
for ti=1:NIter
    x=X(1,ti);
    y=X(2,ti);
    R_OA=sqrt(x^2+y^2);
    vA=X(3:4,ti);
    beta=atan2(y,x);
    beta_dot=[-y, x]*vA/R_OA^2;
    dF_dx=2*n*sign(x)*abs(x)^(2*n-1)/a^(2*n);
    dF_dy=2*n*sign(y)*abs(y)^(2*n-1)/b^(2*n);
    d2F_dxy=0;
    d2F_dxx=2*n*(2*n-1)*abs(x)^(2*n-2)/a^(2*n);
    d2F_dyy=2*n*(2*n-1)*abs(y)^(2*n-2)/b^(2*n);
    kappa(ti)=(-dF_dy^2*d2F_dxx+2*dF_dx*dF_dy*d2F_dxy-dF_dx^2*d2F_dyy)/(dF_dx^2+dF_dy^2)^(3/2);
    if n~=1
        dF_dx_dot=2*n*(2*n-1)*abs(x)^(2*n-2)*vA(1)/a^(2*n);
        dF_dy_dot=2*n*(2*n-1)*abs(y)^(2*n-2)*vA(2)/b^(2*n);
    else
        dF_dx_dot=2*vA(1)/a^(2*n);
        dF_dy_dot=2*vA(2)/b^(2*n);
    end
    beta_bar=atan2(dF_dx,-dF_dy);
    if beta~=0
        %beta_bar_dot=(b^(2*n)*(2*n-1)*sign(cos(beta))*abs(cos(beta))^(2*n-2))/(a^(2*n)*sign(sin(beta))*abs(sin(beta))^(2*n))*beta_dot;
        beta_bar_dot=(dF_dx*dF_dy_dot-dF_dy*dF_dx_dot)/(dF_dx^2+dF_dy^2);
    else
        beta_bar_dot=0;
    end
    f(1:2,1)=X(3:4,ti);
    f(3:4,1)=v0*beta_bar_dot*[-sin(beta_bar);cos(beta_bar)];
    norm_f(ti)=norm(f);
    beta_bar_dot_arr(ti)=beta_bar_dot;
    X(:,ti+1)=X(:,ti)+dt*f;
    time(ti+1)=time(ti)+dt;
end
norm_f(ti+1)=norm_f(ti);
beta_bar_dot_arr(ti+1)=beta_bar_dot_arr(ti);
max_beta_bar_dot=max(beta_bar_dot_arr);
R_circ=v0/max_beta_bar_dot;

figure
hold all
axis equal
rectangle('position',[-w/2,-h/2,w,h])
%rectangule with rounded corners
%rectangle('position',[-w/2-rho_A,-h/2-rho_A,w+2*rho_A,h+2*rho_A])%,'curvature',1/rho_A)
plot(w/2+rho_A*circCos(1:50),h/2++rho_A*circSin(1:50),'k')
plot([w/2,-w/2],rho_A+[h/2,h/2],'k')
plot(-w/2+rho_A*circCos(51:100),h/2++rho_A*circSin(51:100),'k')
plot(-rho_A+[-w/2,-w/2],[-h/2,h/2],'k')
plot(-w/2+rho_A*circCos(101:150),-h/2++rho_A*circSin(101:150),'k')
plot([w/2,-w/2],-rho_A+[-h/2,-h/2],'k')
plot(w/2+rho_A*circCos(151:201),-h/2++rho_A*circSin(151:201),'k')
plot(rho_A+[w/2,w/2],[-h/2,h/2],'k')
%super-ellipse
plot(xsup,ysup,'color',[0.5,0,1],'linewidth',0.5)
%circle
plot(R_circ*circCos,R_circ*circSin,'k--')
%circle
plot(RcO*circCos,RcO*circSin,'g--')
plot(X(1,:),X(2,:),'b--');
hold on;
quiver(X(1,:),X(2,:),X(3,:),X(4,:),'r')


figure
subplot(2,1,1)
plot(time,norm_f)
ylabel('$||f||$')
subplot(2,1,2)
plot(time,beta_bar_dot_arr)
ylabel('$\dot{\bar \beta}$')


%%
% Convex polygon and its centroid
rV=[0,0;5,0;5,5;0,5]';
nV=size(rV,2);
rV=[rV,rV(:,1)]; %Add the first vertex in the list for the calculation purposes

circCos=cos(0:pi/100:2*pi);
circSin=sin(0:pi/100:2*pi);

A=0;
Cx=0;
Cy=0;
for i=1:nV
    temp=rV(1,i)*rV(2,i+1)-rV(2,i)*rV(1,i+1);
    Cx=Cx+(rV(1,i)+rV(1,i+1))*temp;
    Cy=Cy+(rV(2,i)+rV(2,i+1))*temp;
    A=A+temp;
end
A=0.5*A;
Cx=Cx/(6*A);
Cy=Cy/(6*A);

for i=1:nV
    distV(i)=norm(rV(:,i)-[Cx;Cy]);
    mLi=(rV(2,i+1)-rV(2,i))/(rV(1,i+1)-rV(1,i));
    cLi=rV(2,i)-mLi*rV(1,i);
    if isinf(mLi)
        rProj(1,1)=rV(1,i);
        rProj(2,1)=Cy;
    else
        rProj(1,1)=(mLi*Cy+Cx-mLi*cLi)/(1+mLi^2);
        rProj(2,1)=mLi*rProj(1,1)+cLi;
    end
    distL(i)=norm(rProj-[Cx;Cy]);
end

R_inscr=min([distV, distL]);
R_circum=max(distV);
PD_max=R_circum/R_inscr

figure
hold all
plot(rV(1,:),rV(2,:))
plot(Cx,Cy,'ro')
plot(Cx+R_inscr*circCos,Cy+R_inscr*circSin)
plot(Cx+R_circum*circCos,Cy+R_circum*circSin);


%%
%Circles and their tangents (Finding using their parametric reprenstation)
R1=5;
R2=3;
C1=[0,-0]';
C2=[0,20]';
lambda=0:0.005:1;
circCos=cos(0:pi/100:2*pi);
circSin=sin(0:pi/100:2*pi);

for i=1:length(lambda)
    lambda1=lambda(i);
    for j=1:length(lambda)
        lambda2=lambda(j);
        r1=C1+R1*[cos(2*pi*lambda1);sin(2*pi*lambda1)];
        r2=C2+R2*[cos(2*pi*lambda2);sin(2*pi*lambda2)];
        t1=[cos(2*pi*lambda1+pi/2), sin(2*pi*lambda1+pi/2)]';
        t2=[cos(2*pi*lambda2+pi/2), sin(2*pi*lambda2+pi/2)]';
        m1=tan(2*pi*lambda1+pi/2);
        m2=tan(2*pi*lambda2+pi/2);
        theta=atan2(r2(2)-r1(2),r2(1)-r1(1));
        m=tan(theta);
        t=[cos(theta),sin(theta)]';
        dc=cross([t1;0],[t;0])-cross([t2;0],[t;0]);
        G1(i,j)=dc(3);
        dc=cross([t1;0],[t;0])-cross([t2;0],[t;0]);
        G2(i,j)=dc(3);
        G(i,j)=(m-m1)^2+(m-m2)^2;
        G(i,j)=(tanh(m)-tanh(m1))^2+(tanh(m)-tanh(m2))^2;
        Gm=40;
        if G(i,j)>Gm
            G(i,j)=Gm;
        elseif G(i,j)<-Gm
            G(i,j)=-Gm;
        end
        if m>Gm
            M(i,j)=Gm;
        else
            M(i,j)=m;
        end
        if m1>Gm
            M1(i,j)=Gm;
        else
            M1(i,j)=m1;
        end
        if m2>Gm
            M2(i,j)=Gm;
        else
            M2(i,j)=m2;
        end
        
        %G2(i,j)=t1'*t+t2'*t;
    end
end
theta=@(x) atan2(C2(2)+R2*sin(2*pi*x(2))-(C1(2)+R1*sin(2*pi*x(1))),  C2(1)+R2*cos(2*pi*x(2))-(C1(1)+R1*cos(2*pi*x(1))));
m0=@(x) (C2(2)+R2*sin(2*pi*x(2))-(C1(2)+R1*sin(2*pi*x(1))))/(C2(1)+R2*cos(2*pi*x(2))-(C1(1)+R1*cos(2*pi*x(1))));
%dCross= [cos(2*pi*x(1)+pi/2);sin(2*pi*x(1)+pi/2)]'*[cos(theta(x));sin(theta(x))]-[cos(2*pi*x(2)+pi/2);sin(2*pi*x(2)+pi/2)]'*[cos(theta(x));sin(theta(x))];
%f=@(x) [cos(2*pi*x(1)+pi/2)*sin(theta(x))-sin(2*pi*x(1)+pi/2)*cos(theta(x))]+[cos(2*pi*x(2)+pi/2)*sin(theta(x))-sin(2*pi*x(2)+pi/2)*cos(theta(x))];
f=@(x) (m0(x)-tan(2*pi*x(1)+pi/2))^2+ (m0(x)-tan(2*pi*x(2)+pi/2))^2
A=[1,0;0,1;-1,0;0,-1];
b=[1,1,0,0]';
g(:,1)=fmincon(f,[0.7,0.7],A,b)
g(:,2)=fmincon(f,[0.75,0.45],A,b)
g(:,3)=fmincon(f,[0.24,0.7],A,b)
g(:,4)=fmincon(f,[0.24,0.24],A,b)

X=-10:0.1:25;
for i=1:4
    r1(:,i)=C1+R1*[cos(2*pi*g(1,i));sin(2*pi*g(1,i))];
    m1(i)=tan(2*pi*g(1,i)+pi/2);
    c1(i)=r1(2,i)-m1(i)*r1(1,i);
    Y(i,:)=m1(i).*X+c1(i)
end

figure
plot(C1(1)+R1*circCos,C1(2)+R1*circSin,'b')
hold on;
plot(C2(1)+R2*circCos,C2(2)+R2*circSin,'b')
for i=1:4
    plot(X,Y(i,:));
end
% syms g1 g2
% m00=(C2(2)+R2*sin(2*pi*g2)-(C1(2)+R1*sin(2*pi*g1)))/(C2(1)+R2*cos(2*pi*g2)-(C1(1)+R1*cos(2*pi*g1)));
% f=(m00-tan(2*pi*g1+pi/2))^2 + (m00-tan(2*pi*g2+pi/2))^2;
% gradF=gradient(f,[g1,g2]);
% x0=solve(gradF==0)
figure
surface(G)
if  (0)
    figure
    surface(G1)
    
    figure
    surface(G2);
    
    figure
    surface(G)
    
    figure
    surface(M);
    hold on
    surface(M1);
    surface(M2);
end


%% For rectangle with rounded corners and common tangents
rV={[0,1;6,2;6,5;4,9;]',[10,3;14,2;15,7;10.5,8]',[2,-10;7,-12;9,-5;6,-4;]',[12,15;14,18;11,19]'};
%rV={[0,1;6,2;6,5;4,9;2,11]',[1,-10;3,-10;8,-7;5,-4;2,-6]'};
%rV={[10,3;15,2;13,5;]',[0,-10;3,-10;8,-7;5,-4;1,-6]'};
%rV={[0,0;1,0;1,1;0,1]',[5,0;7,0;6,1]'}
%rV={[10,3;15,2;13,5;]',[12,15;14,20;11,19]'};
%rV={[0,-10;10,-10;10,0;0,0;]',[15,3;20,3;20,8;15,8]'};
rho_D=1;
colors={[0,0,1],[1,0,0],[0,0,0],[0.5,0.5,0.5],[0,0.5,1],[1,0.5,1],[0.5,0.5,1],[1,.1,0.5],[0.6,0.4,0.5],[0.1,0.1,0.4],[0,0.5,0.4]};
fontSize=18;

options=optimoptions('fmincon','display','off','algorithm','sqp-legacy');

NO=length(rV);
g=sym('g',[NO,1]);
%syms g
figure
hold all
for k=1:NO
    nV(k)=size(rV{k},2);
    posV=rV{k};
    posV=[posV,posV(:,1)];
    rV{k}=posV;
    plot(posV(1,:),posV(2,:),'color',colors{3})
    
    %Calculate the vertices (outer vertices) of the approximated obstacles
    countV=0;
    for i=1:nV(k)
        mLi=(posV(2,i+1)-posV(2,i))/(posV(1,i+1)-posV(1,i));
        if mLi==0
            xi1=posV(1,i);
            xi2=posV(1,i);
            yi1=posV(2,i)+rho_D;
            yi2=posV(2,i)-rho_D;
        elseif isinf(mLi)
            xi1=posV(1,i)+rho_D;
            xi2=posV(1,i)-rho_D;
            yi1=posV(2,i);
            yi2=posV(2,i);
        else
            dx=sqrt(rho_D^2/(1+1/mLi^2));
            xi1=posV(1,i)+dx;
            xi2=posV(1,i)-dx;
            yi1=posV(2,i)+dx*(-1/mLi);
            yi2=posV(2,i)-dx*(-1/mLi);
        end
        
        crossProd=cross([posV(1:2,i+1)-posV(1:2,i);0],[[xi1;yi1]-posV(1:2,i);0]);
        countV=countV+1;
        if crossProd(3)<0
            posV2(:,countV)=[xi1,yi1]';
            
            countV=countV+1;
            if mLi==0
                xip=posV(1,i+1);
                yip=posV(2,i+1)+rho_D;
            elseif isinf(mLi)
                xip=posV(1,i+1)+rho_D;
                yip=posV(2,i+1);
            else
                xip=posV(1,i+1)+dx;
                yip=posV(2,i+1)+dx*(-1/mLi);
            end
            
        else
            posV2(:,countV)=[xi2,yi2]';
            
            countV=countV+1;
            if mLi==0
                xip=posV(1,i+1);
                yip=posV(2,i+1)-rho_D;
            elseif isinf(mLi)
                xip=posV(1,i+1)-rho_D;
                yip=posV(2,i+1);
            else
                xip=posV(1,i+1)-dx;
                yip=posV(2,i+1)-dx*(-1/mLi);
            end
        end
        %corresponding to the next vertex
        posV2(:,countV)=[xip,yip]';
    end
    %Shift the outer vertex corresponding to the first inner vertex to
    %appropriate position
    posV2=posV2(:,[countV,1:countV-1]);
    rV2{k}=posV2;
    nV2=countV;
    posV2=[posV2,posV2(:,1)];
    %Find the perimeter and Plot the outer vertices
    PeriO(k)=0;
    for ii=1:nV(k)
        %Get the inner vertex coordinates
        xV=posV(1,ii);
        yV=posV(2,ii);
        
        %Find angles corresponding to the outer vertices w.r.t the inner
        %one
        ang1=atan2(posV2(2,2*ii-1)-yV,posV2(1,2*ii-1)-xV);
        if ang1<0
            ang1=ang1+2*pi;
        end
        
        ang2=atan2(posV2(2,2*ii)-yV,posV2(1,2*ii)-xV);
        if ang2<0
            ang2=ang2+2*pi;
        end
        if ang2<ang1
            ang2=ang2+2*pi;
        end
        
        %Find the perimeters of the segments of the approximated obstacle
        PeriSOk{k}(2*ii-1,1)=rho_D*(ang2-ang1);
        PeriSOk{k}(2*ii,1)=norm(posV2(:,2*ii)-posV2(:,2*ii+1));
        
        AngOk{k}(2*ii-1,1)=ang1;
        AngOk{k}(2*ii,1)=ang2;
        %plot the circular arc
        plot(xV+rho_D*cos(ang1:(ang2-ang1)/50:ang2),yV+rho_D*sin(ang1:(ang2-ang1)/50:ang2),'color',colors{4})
        %plot the straight line
        plot(posV2(1,2*ii:2*ii+1),posV2(2,2*ii:2*ii+1),'color',colors{4})%,'color',colors{3},'--')
    end
    %Find the total perimeter of the boundary of the approximated obstacle
    PeriO(k)=sum(PeriSOk{k});
    GammaOk{k}(1,1)=0;
    for ii=1:nV(k)
        GammaOk{k}(2*ii,1)=GammaOk{k}(2*ii-1,1)+PeriSOk{k}(2*ii-1)/PeriO(k);
        GammaOk{k}(2*ii+1,1)=GammaOk{k}(2*ii,1)+PeriSOk{k}(2*ii)/PeriO(k);
        GammaCOk{k}(2*ii-1,1)=AngOk{k}(2*ii-1,1)/2/pi;
        GammaCOk{k}(2*ii,1)=AngOk{k}(2*ii,1)/2/pi;
    end
    
    %gk=G(k);
    % syms xOk(gk) yOk(gk)
    xOk{k}=@(g) 0;
    yOk{k}=@(g) 0;
    xCOk{k}=@(g) 0;
    yCOk{k}=@(g) 0;
    dxOk_dg{k}=@(g) 0;
    dyOk_dg{k}=@(g) 0;
    for ii=1:nV(k)
        i1=2*ii-1;
        i2=2*ii;
        i3=2*ii+1;
        %Get the inner vertex coordinates
        xV=posV(1,ii);
        yV=posV(2,ii);
        gk1=GammaOk{k}(i1);
        gk2=GammaOk{k}(i2);
        gk3=GammaOk{k}(i3);
        ang1=AngOk{k}(2*ii-1,1);
        ang2=AngOk{k}(2*ii,1);
        
        if ii==nV(k)
            xOk{k}=@(g) xOk{k}(g)+(gk1 <= g(k) & g(k) <= gk2)*( xV + rho_D*cos(ang1+(g(k)-gk1)*(ang2-ang1)/(gk2-gk1)))+...
                (gk2 < g(k) & g(k) <= gk3)*( posV2(1,i2)+(g(k)-gk2)*(posV2(1,i3)-posV2(1,i2))/(gk3-gk2));
            
            yOk{k}=@(g) yOk{k}(g)+(gk1 <= g(k) & g(k) <= gk2)*( yV + rho_D*sin(ang1+(g(k)-gk1)*(ang2-ang1)/(gk2-gk1)))+...
                (gk2 < g(k) & g(k) <= gk3)*( posV2(2,i2)+(g(k)-gk2)*(posV2(2,i3)-posV2(2,i2))/(gk3-gk2));
            
            xCOk{k}=@(g) xCOk{k}(g)+( GammaCOk{k}(i1) <= g(k) & g(k) <= GammaCOk{k}(i2))*( xV + rho_D*cos(2*pi*g(k)));
            
            yCOk{k}=@(g) yCOk{k}(g)+(GammaCOk{k}(i1) <= g(k) & g(k) <= GammaCOk{k}(i2))*( yV + rho_D*sin(2*pi*g(k)));
            
            dxOk_dg{k}=@(g) dxOk_dg{k}(g)+(gk1 <= g(k) & g(k) <= gk2)*(-rho_D*sin(ang1+(g(k)-gk1)*(ang2-ang1)/(gk2-gk1)))*(ang2-ang1)/(gk2-gk1)+...
                (gk2 < g(k) & g(k) <= gk3)*((posV2(1,i3)-posV2(1,i2))/(gk3-gk2));
            
            dyOk_dg{k}=@(g) dyOk_dg{k}(g)+(gk1 <= g(k) & g(k) <= gk2)*(rho_D*cos(ang1+(g(k)-gk1)*(ang2-ang1)/(gk2-gk1)))*(ang2-ang1)/(gk2-gk1)+...
                (gk2 < g(k) & g(k) <= gk3)*((posV2(2,i3)-posV2(2,i2))/(gk3-gk2));
        else
            xOk{k}=@(g) xOk{k}(g)+(gk1 <= g(k) & g(k) <= gk2)*( xV + rho_D*cos(ang1+(g(k)-gk1)*(ang2-ang1)/(gk2-gk1)))+...
                (gk2 < g(k) & g(k) < gk3)*( posV2(1,i2)+(g(k)-gk2)*(posV2(1,i3)-posV2(1,i2))/(gk3-gk2));
            
            yOk{k}=@(g) yOk{k}(g)+(gk1 <= g(k) & g(k) <= gk2)*( yV + rho_D*sin(ang1+(g(k)-gk1)*(ang2-ang1)/(gk2-gk1)))+...
                (gk2 < g(k) & g(k) < gk3)*( posV2(2,i2)+(g(k)-gk2)*(posV2(2,i3)-posV2(2,i2))/(gk3-gk2));
            
            xCOk{k}=@(g) xCOk{k}(g)+( GammaCOk{k}(i1) <= g(k) & g(k) < GammaCOk{k}(i2))*( xV + rho_D*cos(2*pi*g(k)));
            
            yCOk{k}=@(g) yCOk{k}(g)+(GammaCOk{k}(i1) <= g(k) & g(k) < GammaCOk{k}(i2))*( yV + rho_D*sin(2*pi*g(k)));
            
            dxOk_dg{k}=@(g) dxOk_dg{k}(g)+(gk1 <= g(k) & g(k) <= gk2)*(-rho_D*sin(ang1+(g(k)-gk1)*(ang2-ang1)/(gk2-gk1)))*(ang2-ang1)/(gk2-gk1)+...
                (gk2 < g(k) & g(k) < gk3)*((posV2(1,i3)-posV2(1,i2))/(gk3-gk2));
            
            dyOk_dg{k}=@(g) dyOk_dg{k}(g)+(gk1 <= g(k) & g(k) <= gk2)*(rho_D*cos(ang1+(g(k)-gk1)*(ang2-ang1)/(gk2-gk1)))*(ang2-ang1)/(gk2-gk1)+...
                (gk2 < g(k) & g(k) < gk3)*((posV2(2,i3)-posV2(2,i2))/(gk3-gk2));
        end
        
        Px{k}(ii)=(gk1 <= g(k) & g(k) < gk2)*( xV + rho_D*cos(ang1+(g(k)-gk1)*(ang2-ang1)/(gk2-gk1)))+...
            (gk2 <= g(k) & g(k) < gk3)*( posV2(1,i2)+(g(k)-gk2)*(posV2(1,i3)-posV2(1,i2))/(gk3-gk2));
        dPx_dg{k}(ii)=(gk1 <= g(k) & g(k) < gk2)*(-rho_D*sin(ang1+(g(k)-gk1)*(ang2-ang1)/(gk2-gk1)))*(ang2-ang1)/(gk2-gk1)+...
            (gk2 <= g(k) & g(k) < gk3)*((posV2(1,i3)-posV2(1,i2))/(gk3-gk2));
        
        d2Px_dg2{k}(ii)=(gk1 <= g(k) & g(k) < gk2)*(-rho_D*cos(ang1+(g(k)-gk1)*(ang2-ang1)/(gk2-gk1)))*((ang2-ang1)/(gk2-gk1))^2;
        
        Py{k}(ii)=(gk1 <= g(k) & g(k) < gk2)*( yV + rho_D*sin(ang1+(g(k)-gk1)*(ang2-ang1)/(gk2-gk1)))+...
            (gk2 <= g(k) & g(k) < gk3)*( posV2(2,i2)+(g(k)-gk2)*(posV2(2,i3)-posV2(2,i2))/(gk3-gk2));
        
        dPy_dg{k}(ii)=(gk1 <= g(k) & g(k) < gk2)*(rho_D*cos(ang1+(g(k)-gk1)*(ang2-ang1)/(gk2-gk1)))*(ang2-ang1)/(gk2-gk1)+...
            (gk2 <= g(k) & g(k) < gk3)*((posV2(2,i3)-posV2(2,i2))/(gk3-gk2));
        %xOk=matlabfunction(xOk,'Vars',gk);
        d2Py_dg2{k}(ii)=(gk1 <= g(k) & g(k) < gk2)*(-rho_D*sin(ang1+(g(k)-gk1)*(ang2-ang1)/(gk2-gk1)))*((ang2-ang1)/(gk2-gk1))^2;
        %         funPx{k}(1,4*ii-3)=gk1 <= gk <= gk2;
        %         funPx{k}(1,4*ii-2)= xV + rho_D*cos(ang1+(gk-gk1)*(ang2-ang1)/(gk2-gk1));
        %         funPx{k}(1,4*ii-1)=gk2 <= gk <= gk3;
        %         funPx{k}(1,4*ii)=posV2(1,i2)+(gk-gk2)*(posV2(1,i3)-posV2(1,i2))/(gk3-gk2);
        %
        %         funPy{k}(1,4*ii-3)=gk1 <= gk <= gk2;
        %         funPy{k}(1,4*ii-2)= yV + rho_D*sin(ang1+(gk-gk1)*(ang2-ang1)/(gk2-gk1));
        %         funPy{k}(1,4*ii-1)=gk2 <= gk <= gk3;
        %         funPy{k}(1,4*ii)=posV2(2,i2)+(gk-gk2)*(posV2(2,i3)-posV2(2,i2))/(gk3-gk2);
    end
    %    xOk{k}=sum(Px{k});
    %    yOk{k}=sum(Py{k});
    %    dxOk_dg{k}=sum(dPx_dg{k});
    %    dyOk_dg{k}=sum(dPy_dg{k});
    %    d2xOk_dg2{k}=sum(d2Px_dg2{k});
    %    d2yOk_dg2{k}=sum(d2Py_dg2{k});
    
    %     xOk=matlabFunction(piecewise(funPx{k}),'Vars',[gk]);
    %     yOk=matlabFunction(piecewise(funPy{k}),'Vars',[gk]);
    
    clear posV2;
    clear posV;
end

A=[eye(NO);-eye(NO)];
%b=[1,1,0,0]';

gamma=0:.01:1;
countP=0;
countTotT=0;
countTO=zeros(NO,1);
for k=1:NO-1
    for kk=k+1:NO
        countP=countP+1;
        m=@(g) (yOk{k}(g) - yOk{kk}(g))/(xOk{k}(g) - xOk{kk}(g));
        m1=@(g) dyOk_dg{k}(g)/dxOk_dg{k}(g);
        m2=@(g) dyOk_dg{kk}(g)/dxOk_dg{kk}(g);
        f=@(g) (m(g)-m1(g))^2+(m(g)-m2(g))^2;
        
        mC=@(g)  (yCOk{k}(g) - yCOk{kk}(g))/(xCOk{k}(g) - xCOk{kk}(g));
        mC1=@(g) tan(2*pi*g(k)+pi/2);
        mC2=@(g) tan(2*pi*g(kk)+pi/2);
        fC=@(g) (mC(g)-mC1(g))^2+(mC(g)-mC2(g))^2;
        
        countT=0;
        
        for i=1:nV(k)
            for j=1:nV(kk)
                b=zeros(2*NO,1);
                b([k,kk],1)=[GammaOk{k}(2*i), GammaOk{kk}(2*j)]';
                b(NO+[k,kk],1)=[-GammaOk{k}(2*i-1), -GammaOk{kk}(2*j-1)]';
                g0=zeros(1,NO);
                alp1=.95;
                g0(k)=alp1*GammaOk{k}(2*i-1)+(1-alp1)*GammaOk{k}(2*i);
                g0(kk)=alp1*GammaOk{kk}(2*j-1)+(1-alp1)*GammaOk{kk}(2*j);
                gt0=fmincon(f,g0,A,b,[],[],[],[],[],options);
                ftol=1e-7;
                gtol=1e-5;
                % gt0=fmincon(f,[GammaOk{k}(2*i-1),GammaOk{kk}(2*j-1)],A,b);
                if f(gt0)<ftol && GammaOk{k}(2*i-1)<gt0(k) && gt0(k)<GammaOk{k}(2*i) && GammaOk{kk}(2*j-1)<gt0(kk) && gt0(kk)<GammaOk{kk}(2*j)
                    if countT==0
                        countT=countT+1;
                        gt(:,countT)=gt0;
                    elseif ~sum( (abs(gt(k,:)-gt0(k))<gtol) & (abs(gt(kk,:)-gt0(kk))<gtol))
                        countT=countT+1;
                        gt(:,countT)=gt0;
                    end
                end
                alp2=0.24;
                g0(k)=alp2*GammaOk{k}(2*i-1)+(1-alp2)*GammaOk{k}(2*i);
                g0(kk)=alp2*GammaOk{kk}(2*j-1)+(1-alp2)*GammaOk{kk}(2*j);
                gt0=fmincon(f,g0,A,b,[],[],[],[],[],options);
                % gt0=fmincon(f,[GammaOk{k}(2*i-1),GammaOk{kk}(2*j-1)],A,b);
                if f(gt0)<ftol && GammaOk{k}(2*i-1)<gt0(k) && gt0(k)<GammaOk{k}(2*i) && GammaOk{kk}(2*j-1)<gt0(kk) && gt0(kk)<GammaOk{kk}(2*j)
                    if countT==0
                        countT=countT+1;
                        gt(:,countT)=gt0;
                    elseif ~sum( (abs(gt(k,:)-gt0(k))<gtol) & (abs(gt(kk,:)-gt0(kk))<gtol))
                        countT=countT+1;
                        gt(:,countT)=gt0;
                    end
                end
                
                g0(k)=alp1*GammaOk{k}(2*i-1)+(1-alp1)*GammaOk{k}(2*i);
                g0(kk)=alp2*GammaOk{kk}(2*j-1)+(1-alp2)*GammaOk{kk}(2*j);
                gt0=fmincon(f,g0,A,b,[],[],[],[],[],options);
                % gt0=fmincon(f,[GammaOk{k}(2*i-1),GammaOk{kk}(2*j-1)],A,b);
                if f(gt0)<ftol && GammaOk{k}(2*i-1)<gt0(k) && gt0(k)<GammaOk{k}(2*i) && GammaOk{kk}(2*j-1)<gt0(kk) && gt0(kk)<GammaOk{kk}(2*j)
                    if countT==0
                        countT=countT+1;
                        gt(:,countT)=gt0;
                    elseif ~sum( (abs(gt(k,:)-gt0(k))<gtol) & (abs(gt(kk,:)-gt0(kk))<gtol))
                        countT=countT+1;
                        gt(:,countT)=gt0;
                    end
                end
                
                g0(k)=alp2*GammaOk{k}(2*i-1)+(1-alp2)*GammaOk{k}(2*i);
                g0(kk)=alp1*GammaOk{kk}(2*j-1)+(1-alp1)*GammaOk{kk}(2*j);
                gt0=fmincon(f,g0,A,b,[],[],[],[],[],options);
                % gt0=fmincon(f,[GammaOk{k}(2*i-1),GammaOk{kk}(2*j-1)],A,b);
                if f(gt0)<ftol && GammaOk{k}(2*i-1)<gt0(k) && gt0(k)<GammaOk{k}(2*i) && GammaOk{kk}(2*j-1)<gt0(kk) && gt0(kk)<GammaOk{kk}(2*j)
                    if countT==0
                        countT=countT+1;
                        gt(:,countT)=gt0;
                    elseif ~sum( (abs(gt(k,:)-gt0(k))<gtol) & (abs(gt(kk,:)-gt0(kk))<gtol))
                        countT=countT+1;
                        gt(:,countT)=gt0;
                    end
                end
                
            end
        end
        
        
        %         dm1_dgk=d2yOk_dg2{k}/dxOk_dg{k}-dyOk_dg{k}*d2xOk_dg2{k}/dxOk_dg{k}^2;
        %         dm2_dgkk=d2yOk_dg2{kk}/dxOk_dg{kk}-dyOk_dg{kk}*d2xOk_dg2{kk}/dxOk_dg{kk}^2;
        %         dm_dgk=dyOk_dg{k}/(xOk{k}-xOk{kk}) - (yOk{k}-yOk{kk})*dxOk_dg{k}/(xOk{k}-xOk{kk})^2;
        %         dm_dgkk=-dyOk_dg{kk}/(xOk{k}-xOk{kk}) + (yOk{k}-yOk{kk})*dxOk_dg{kk}/(xOk{k}-xOk{kk})^2;
        
        %         df_dgk=2*(m-m1)*(dm_dgk-dm1_dgk)+2*(m-m2)*dm_dgk;
        %         df_dgkk=2*(m-m1)*(dm_dgkk)+2*(m-m2)*(dm_dgkk-dm2_dgkk);
        
        %f=matlabFunction(f);
        %         df_dgk=matlabFunction(df_dgk);
        %         df_dgkk=matlabFunction(df_dgkk);
        if(0)
            for i=1:length(gamma)
                x1(i)=xOk{k}([gamma(i),0]);
                y1(i)=yOk{k}([gamma(i),0]);
                for j=1:length(gamma)
                    x2(j)=xOk{kk}([0,gamma(j)]);
                    y2(j)=yOk{kk}([0,gamma(j)]);
                    
                    F(i,j)=f([gamma(i),gamma(j)]);
                    M(i,j)=m([gamma(i),gamma(j)]);
                    M1(i,j)=m1([gamma(i),gamma(j)]);
                    M2(i,j)=m2([gamma(i),gamma(j)]);
                    
                    if F(i,j)>20
                        F(i,j)=20;
                    elseif F(i,j)<-20
                        F(i,j)=-20;
                    end
                    Gm=500;
                    if M1(i,j)>Gm
                        M1(i,j)=Gm;
                    elseif M1(i,j)<-Gm
                        M1(i,j)=-Gm;
                    end
                    
                    if M2(i,j)>Gm
                        M2(i,j)=Gm;
                    elseif M2(i,j)<-Gm
                        M2(i,j)=-Gm;
                    end
                    %dF_dk(i,j)=df_dgk
                end
            end
        end
        
        
        %X=-2:0.1:15;
        for i=1:size(gt,2)
            r1(:,i)=[xOk{k}([gt(:,i)']),yOk{k}(gt(:,i)')];
            r2(:,i)=[xOk{kk}([gt(:,i)']),yOk{kk}(gt(:,i)')];
            
            % Store the tangent points on the curves
            countTO(k)= countTO(k)+1;
            countTO(kk)= countTO(kk)+1;
            rTO{k}(:,countTO(k))=r1(:,i);
            rTO{kk}(:,countTO(kk))=r2(:,i);
            gTO{k}(countTO(k))=gt(k,i);
            gTO{kk}(countTO(kk))=gt(kk,i);
        end
        hold on;
        for i=1:size(gt,2)
            plot([r1(1,i),r2(1,i)],[r1(2,i),r2(2,i)],'color',colors{countP+4});
        end
        
        
        
        clear gt
        countTotT=countTotT+countT;
    end
    
end




% gt(:,1)=fmincon(f,[0.7,0.7],A,b)
% gt(:,2)=fmincon(f,[0.6,0.21],A,b)
% gt(:,3)=fmincon(f,[0.24,0.6],A,b)
% gt(:,4)=fmincon(f,[0.24,0.24],A,b)

% gt(:,1)=fsolve(f,[0.7,0.7])
% gt(:,2)=fsolve(f,[0.55,0.1])
% gt(:,3)=fsolve(f,[0.34,0.05])
% gt(:,4)=fsolve(f,[0.14,0.24])


figure
surface(F)

figure
surface(M1)


%%
% f for a specific range of gamma
clear F
gamma1=GammaOk{1}(1):(GammaOk{1}(1)+GammaOk{1}(2))/300:GammaOk{1}(2);
gamma2=GammaOk{2}(1):(GammaOk{2}(1)+GammaOk{2}(2))/300:GammaOk{2}(2);

for i=1:length(gamma1)
    x1(i)=xOk{k}([gamma1(i),0]);
    y1(i)=yOk{k}([gamma1(i),0]);
    for j=1:length(gamma2)
        x2(j)=xOk{kk}([0,gamma2(j)]);
        y2(j)=yOk{kk}([0,gamma2(j)]);
        
        F(i,j)=f([gamma1(i),gamma2(j)]);
        M(i,j)=m([gamma1(i),gamma2(j)]);
        M1(i,j)=m1([gamma1(i),gamma2(j)]);
        M2(i,j)=m2([gamma1(i),gamma2(j)]);
        if F(i,j)>1200
            F(i,j)=1200;
        elseif F(i,j)<-1200
            F(i,j)=-1200;
        end
        if M1(i,j)>Gm
            M1(i,j)=Gm;
        elseif M1(i,j)<-Gm
            M1(i,j)=-Gm;
        end
        
        if M2(i,j)>Gm
            M2(i,j)=Gm;
        elseif M2(i,j)<-Gm
            M2(i,j)=-Gm;
        end
        %dF_dk(i,j)=df_dgk
    end
end

figure
surface(F)



%%
r1=[1,-10]';
r2=[10,15]';
dr=r2-r1;
drx=r2(1)-r1(1);
dry=r2(2)-r1(2);
mL1=(r2(2)-r1(2))/(r2(1)-r1(1));
cL1=r1(2)-mL1*r1(1);
dr_hat=dr/norm(dr);
%check with the lines between the inner vertices
for i=1:nVO(k)
    mL2=(rVO{k}(2,i)-rVO{k}(2,i+1))/(rVO{k}(1,i)-rVO{k}(1,i+1));
    cL2=rVO{k}(2,i)-mL2*rVO{k}(1,i);
    drx2=rVO{k}(1,i+1)-rVO{k}(1,i);
    x_int=(cL2-cL1)/(mL1-mL2);
    %y_int=mL1*x_int+cL1;
    lambda1=(x_int-r1(1))/drx;
    lambda2=(x_int-rVO{k}(1,i))/drx2;
    if (0<= lambda1) && (lambda1<=1) && (0<= lambda2) && (lambda2<=1) % %if line segments intersects the line joining two vertices
        countRemoveTang=countRemoveTang+1;
        indRemoveTang(countRemoveTang)=kt;
        flagRemoveTang=1;
        break;
    end
    
end

g00=0:0.01:1;
g0=zeros(1,NO);
for j=1:length(gamma)
    g0(k)=g00(j);
    
    F(j)=f(g0)^2;
end


%%
%Test double integrator dynamics constrained on a circular path of radius R
R=10;
X(:,1)=[10,0]';
dt=0.005;
t(1)=0;
s(1)=0;
s_dot(1)=0;
NSim=10000;
for i=1:NSim
    qx=.01;
    qy=.01;
    s_ddot(i)=((qx+qy)+1/R*(cos(s(i)/R)+sin(s(i)/R))*s_dot(i)^2)/(cos(s(i)/R)-sin(s_dot(i)/R));
    %     s_ddot=(qx+1/R*(cos(s(i)/R)*s_dot(i)^2)/(-sin(s(i)/R)));
    %     s_ddot=sqrt(s_dot(i)^4/R^2-qx^2-qy^2);
    s_dot(i+1)=s_dot(i)+dt*s_ddot(i);
    s(i+1)=s(i)+dt*s_dot(i);
    X(:,i+1)=R*[cos(s(i+1)/R),sin(s(i+1)/R)]';
    X_ddot(:,i+1)=[-sin(s(i+1)/R),cos(s(i+1)/R)]'*s_ddot(i)-1/R*[cos(s(i+1)/R),sin(s(i+1)/R)]'*s_dot(i+1)^2;
end

figure
plot(X(1,:),X(2,:))

figure
plot(t,X_ddot(1,:),t,X_ddot(2,:))


%%
%Test constraints of speeds on a path with straight-circ-straight
R=15;
qm=1;
vm=4;
vmc=sqrt(qm*R);
s1=.1*R; s3=.1*R;  %less than R/2
s2=.1*R;
v1m=sqrt(2*qm*s1);
v2m=sqrt(2*qm*s3);
v1=0:v1m/100:v1m;
v2=0:v2m/100:v2m;

for i=1:length(v1)
    v21(i)=sqrt(v1(i)^2+2*qm*max(v1(i)))/(1+2*s2/R);
    v22(i)=sqrt(v1(i)^2-2*s2*(qm-v1(i)^2/R));
    v21(i)=sqrt(v1(i)^2+2*s2*(qm));
    v22(i)=sqrt(2*s2*(qm)-v1(i)^2);
end

figure
rectangle('position',[0,0,v1m,v2m])
hold on
plot(v1,v1,'b')
plot(v1,v21,v1,v22);



%Minimum time for straight line 
v1=0:vm/100:vm;
v2=0:vm/100:vm;
for i=1:length(v1)
    for j=1:length(v2)
        if 2*qm*s1>abs(v1(i)^2-v2(j)^2)
            if s1<0.5*(2*vm^2-v1(i)^2-v2(j)^2)/qm
                tau(i,j)=(0.5*(2*(v1(i)^2+v2(j)^2+2*s1*qm))-v1(i))/qm+(0.5*(2*(v1(i)^2+v2(j)^2+2*s1*qm))-v2(j))/qm;
            else
                tau(i,j)=s1/vm-0.5*(2*vm^2-v1(i)^2-v2(j)^2)/(qm*vm)+(vm-v1(i))/qm+(vm-v2(j))/qm;
            end
        else
            tau(i,j)=0;
        end
    end
end

Rb=sqrt(2*vm^2-2*s1*qm);
figure
plot3(Rb*cos(0:pi/100:pi/2),Rb*sin(0:pi/100:pi/2),ones(51,1)*(max(max(tau))))
hold on
surface(v1,v2,tau);

% 3 segment path
vmc99=.99*vmc;
v1=0:vmc99/60:vmc99;
v2=0:vmc99/60:vmc99;
LargeNumber=130;
k2=sqrt(qm/R);
for i=1:length(v1)
    for j=1:length(v2)
        %first segment
        if (1)% 2*qm*s1>abs(v1(i)^2)
            if s1<0.5*(2*vm^2-v1(i)^2)/qm
                t1(i,j)=(0.5*sqrt(2*(v1(i)^2+0+2*s1*qm))-v1(i))/qm+(0.5*sqrt(2*(v1(i)^2+0+2*s1*qm))-0)/qm;
            else
                t1(i,j)=s1/vm-0.5*(2*vm^2-v1(i)^2)/(qm*vm)+(vm-v1(i))/qm+(vm)/qm;
            end
        else
            t1(i,j)=LargeNumber;
        end
        
        %third segment
        if (1) %2*qm*s1>abs(v2(j)^2)
            if s1<0.5*(2*vm^2-v2(j)^2)/qm
                t3(i,j)=(0.5*sqrt(2*(v2(j)^2+0+2*s1*qm))-v2(j))/qm+(0.5*sqrt(2*(v2(j)^2+0+2*s1*qm))-0)/qm;
            else
                t3(i,j)=s1/vm-0.5*(2*vm^2-v2(j)^2)/(qm*vm)+(vm-v2(j))/qm+(vm)/qm;
            end
        else
            t3(i,j)=LargeNumber;
        end
        
        %second segment
        maxV=max(v1(i),v2(j));
        minV=min(v1(i),v2(j));
        if (1)%s2>R/2*(log(vmc^2-minV^2)-log(vmc^2-maxV^2))%R*(log(cosh(atanh(maxV/vmc)))-log(cosh(atanh(minV/vmc))))  %Both expressions are the same
            
            %t_half=1/k2*(acosh(cosh(atanh(v1(i)/vmc))*exp(s2/(2*R)))-atanh(v1(i)/vmc));
            %v_half=vmc*tanh(k2*t_half+atanh(v1(i)/vmc));
            v_half=sqrt(vmc^2-exp(-s2/R)*(vmc^2-minV^2));
            v_half=sqrt(vmc^2-exp(-s2/R)*sqrt((vmc^2-minV^2)*(vmc^2-maxV^2)));
            t_half=1/k2*(atanh(v_half/vmc)-atanh(minV/vmc));
            t_half2=1/k2*(atanh(v_half/vmc)-atanh(maxV/vmc));
            se=R*(log(cosh(k2*t_half2+atanh(v_half/vmc)))-log(cosh(atanh(v_half/vmc))))
            t2(i,j)=t_half+t_half2;
        else
            %             v_half=sqrt(vmc^2-exp(-s2/R)*sqrt((vmc^2-minV^2)*(vmc^2-maxV^2)));
            %             t_half=1/k2*(atanh(v_half/vmc)-atanh(minV/vmc));
            %             t_half2=1/k2*(atanh(v_half/vmc)-atanh(maxV/vmc));
            %             se=R*(log(cosh(k2*t_half2+atanh(v_half/vmc)))-log(cosh(atanh(v_half/vmc))))
            %             t2(i,j)=t_half+t_half2;
            t2(i,j)=LargeNumber;
        end
        tau1(i,j)=t1(i,j)+t2(i,j)+t3(i,j);
    end
end

Rb=sqrt(2*vm^2-2*s1*qm);
figure
%plot3(Rb*cos(0:pi/100:pi/2),Rb*sin(0:pi/100:pi/2),ones(51,1)*(max(max(tau))))
hold on
surface(v1,v2,tau1);

%Rb=sqrt(2*vmc^2-2*s1*qm);
figure
% plot3(Rb*cos(0:pi/100:pi/2),Rb*sin(0:pi/100:pi/2),ones(51,1)*(max(max(tau))))
% hold on
surface(v1,v2,t2);
%%
%hypergeom1=@(y) hypergeom(1/4,1/2,5/4,y);
p=@(x) tan(s2/R+0.5*atan(x(1)^2/sqrt(vmc^4-x(1)^4))+0.5*atan(x(2)^2/sqrt(vmc^4-x(2)^4)));
v_half=@(x) vmc*(p(x)^2/(1+p(x)^2))^(1/4);
f1=@(x) v_half(x)*hypergeom([1/4,1/2],5/4,v_half(x)^4/vmc^4);
f2=@(x) x(1)*hypergeom([1/4,1/2],5/4,x(1)^4/vmc^4);
f3=@(x) x(2)*hypergeom([1/4,1/2],5/4,x(2)^4/vmc^4);
f=@(x) 2*f1(x)-f2(x)-f3(x);
%f=@(x) -x(1)*hypergeom([1/4,1/2],5/4,x(1)^4/vmc^4)-x(2)*hypergeom([1/4,1/2],5/4,x(2)^4/vmc^4)
%f=@(x) (tan(1+0.5*atan(x(1)^2/sqrt(vmc^4-x(1)^4)+0.5*atan(x(2)^2/sqrt(vmc^4-x(2)^4))))^2/(1+tan(1+0.5*atan(x(1)^2/sqrt(vmc^4-x(1)^4)+0.5*atan(x(2)^2/sqrt(vmc^4-x(2)^4))))^2))^(0.25);
 for i=1:length(v1)
    for j=1:length(v2)
        vh(i,j)=v_half([v1(i),v2(j)]);
        mv(i,j)=max(v1(i),v2(j));
        F1(i,j)=f1([v1(i),v2(j)]);
        F2(i,j)=f2([v1(i),v2(j)]);
        F3(i,j)=f3([v1(i),v2(j)]);
        t22(i,j)=2*F1(i,j)-F2(i,j)-F3(i,j);
    end
 end
 figure
 surf(v1,v2,vh)
 hold on
 surf(v1,v2,mv)
 figure
 surf(v1,v2,t22)
 
 figure
 surf(v1,v2,F1);
 hold on
 surf(v1,v2,F2)
 surf(v1,v2,F3)
 
 figure
 surf(v1,v2,2*F1)
 hold on;
 surf(v1,v2,F2+F3)
 
 %%
 %s2=
 v_half=@(x) vmc*sqrt(tanh(s2/R+0.5*atanh(x(1)^2/vmc^2)+0.5*atanh(x(2)^2/vmc^2)))
f=@(x) 1/(2*k2)*(log((v_half(x)+vmc)/(vmc-v_half(x)))+2*atan(v_half(x)/vmc))-1/(4*k2)*(log((x(1)+vmc)/(vmc-x(1)))+2*atan(x(1)/vmc))-1/(4*k2)*(log((x(2)+vmc)/(vmc-x(2)))+2*atan(x(2)/vmc))
%f=@(x) (tan(1+0.5*atan(x(1)^2/sqrt(vmc^4-x(1)^4)+0.5*atan(x(2)^2/sqrt(vmc^4-x(2)^4))))^2/(1+tan(1+0.5*atan(x(1)^2/sqrt(vmc^4-x(1)^4)+0.5*atan(x(2)^2/sqrt(vmc^4-x(2)^4))))^2))^(0.25);
 for i=40%1:length(v1)
    for j=1:length(v2)
        t22(i,j)=f([v1(i),v2(j)]);
    end
 end
 figure
 plot(v2,t22(i,:))
 %surf(v1,v2,t22)


 %%
 %plot qm^2-v^4/R^2
 qm=10;
 R=16;
 vmc=sqrt(qm*R);
 figure
 u1=@(v) sqrt(qm^2-v^4/R^2)
 fplot(u1,[0,vmc])
 hold on;
 u2=@(v) log(qm^2-v^4/R^2+1)*qm/log(qm^2+1)%qm-v^2/R;
 fplot(u2,[0,vmc])
 u3=@(v) qm-v^4/(vmc^2*R);
 fplot(u3,[0,vmc])
 
 figure
 u4=@(v) u1(v)-u2(v)
 u5=@(v) u1(v)-u3(v)
 fplot(u4,[0,vmc])
 hold on 
 fplot(u5,[0,vmc]) 
 
 %% 
 %%Plot v vs s
 SetPlotDefaults;
 
 qm=10;
 R=4;
 vmc=sqrt(qm*R);
 k2=sqrt(qm/R);
 vs=0.5;
 ve=2.5;
 s2min=R/2*(atanh(ve^2/vmc^2)-atanh(vs^2/vmc^2));
 s2=10;  %choose 2*pi*R>s2>R/2*(atanh(ve^2/vmc^2)-atanh(vs^2/vmc^2))
 v_half=vmc*sqrt(tanh(s2/R+0.5*atanh(vs^2/vmc^2)+0.5*atanh(ve^2/vmc^2)))
 t_opt=1/(2*k2)*(log((v_half+vmc)/(vmc-v_half))+2*atan(v_half/vmc))-1/(4*k2)*(log((vs+vmc)/(vmc-vs))+2*atan(vs/vmc))-1/(4*k2)*(log((ve+vmc)/(vmc-ve))+2*atan(ve/vmc))
 dt=0.01;
 v(1)=vs;
 s(1)=0;
 time(1)=0;
NSim=ceil(t_opt/dt);
flag=0;
for i=1:NSim
    if v(i)<=v_half && flag==0
    v_dot(i)=qm-v(i)^4/vmc^2/R;
    t_half=time(i);
    else
        flag=1;        
        v_dot(i)=-(qm-v(i)^4/vmc^2/R);
    end
    v(i+1)=v(i)+dt*v_dot(i);
    s(i+1)=s(i)+v(i)*dt;
    time(i+1)=time(i)+dt;
end
figure
yyaxis left
plot(time,v)
hold on
plot(time,ones(length(time))*vmc,'.')
plot(time,ones(length(time))*v_half,'--')
plot([t_opt,t_opt],[0.95*ve,0],'--')
text(t_opt,1.05*vmc,'$v_m^c$','fontsize',18)
text(t_half,0.95*v_half,'$\bar{v}$','fontsize',18)
text(0.01,1*vs,'$v_1$','fontsize',18)
text(t_opt,1*ve,'$v_2$','fontsize',18)
text(0.99*t_opt,-0.2,'$t_{opt}$','fontsize',18)
ylabel('$v \;[m/s]$')
hold on
yyaxis right
plot(time,s)
ylabel('$s \;[m]$')
xlabel('$time \;[s]$')

%axis equal

%%
%Find hessian of optimal t for v_dot=qm-v^4/(R^2qm) 
syms v1 v2 S R vm qm
vbar=vm*sqrt(tanh(S/R+0.5*atanh(v1^2/vm^2)+0.5*atanh(v2^2/vm^2)));
t=0.25*sqrt(R/qm)*(2*(log((vm+vbar)/(vm-vbar))+2*atan(vbar/vm))-(log((vm+v1)/(vm-v1))+2*atan(v1/vm))-(log((vm+v2)/(vm-v2))+2*atan(v2/vm)));
nabla2_t=hessian(t,[v1,v2]);
nabla2_t_fun=matlabFunction(nabla2_t, 'vars',[v1,v2,S,R,vm,qm]);
S=1000.1;
qm=10;
R=4;
vm=sqrt(qm*R);
vmc99=.99*vm;
v1=0:vmc99/100:vmc99;
v2=0:vmc99/100:vmc99;
 for i=1:length(v1)
    for j=1:length(v2)
        %vec=mat2cell([v1(i),v2(j),S,R,vm,qm]);
        nabla2_t2{i,j}=nabla2_t_fun(v1(i),v2(j),S,R,vm,qm);
        det_nabla2_t(i,j)=det(nabla2_t2{i,j});
    end
 end
 figure
 surf(v1,v2,det_nabla2_t)
 
 
 %%
 % Plot time for fixed initial speed
 
 
 
 %%
 %Symbolic hessian of tau (actual bounds)
 syms x y s2 R vmc
 p= tan(s2/R+0.5*atan(x^2/sqrt(vmc^4-x^4))+0.5*atan(y^2/sqrt(vmc^4-y^4) ));
v_half=vmc*(p^2/(1+p^2))^(1/4);
f1=v_half*hypergeom([1/4,1/2],5/4,v_half^4/vmc^4);
f2=x*hypergeom([1/4,1/2],5/4,x^4/vmc^4);
f3=y*hypergeom([1/4,1/2],5/4,y^4/vmc^4);
f=matlabFunction(2*f1-f2-f3,'Vars',[x,y,s2,R,vmc]);
