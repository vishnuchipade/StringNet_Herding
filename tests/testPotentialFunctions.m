

%Test Potential func1tion and gradient
% ------------------------------
% Author: Vishnu S. Chipade
% Email: vishnuc@umich.edu
% ------------------------------


alpha=0.5;
c1=0.91;
r0=10;
rm=1;
c0=1;
R_m2_DD=10;
Rjj0=2;
x=0:0.1:R_m2_DD;

% %V(x)=log((c1*x-c2)/(x-xm)+(x-xm)/(c1*x-c2));
c1=(Rjj0-c0)/(R_m2_DD-Rjj0)
c2=R_m2_DD*c1;      
temp1=(c1^2-1)*x.^2-2*c1*c2*x+c2^2-c0^2+2*c0*x;
temp2=(c1^2+1)*x.^2-2*c1*c2*x+c2^2+c0^2-2*c0*x;

V=log((c1*x-c2)./(x-c0)+(x-c0)./(c1*x-c2));
norm_nabla_V=1./x.*(-(c2-c1.*c0)./((c0-x).*(c1.*x-c2)).*temp1./temp2);
x_star1=(c2+c0)/(c1+1)
x_star2=(c2-c0)/(c1-1)

figure
plot(x,V,'b-',x,norm_nabla_V,'b--')
hold on;

%V1(x)=log(k/(x-x0)+(x-x0)/k) 
x=0:0.1:50;
V1=log(Rjj0./(x-c0)+(x-c0)./Rjj0);
norm_nabla_V1=1./x./abs(x-c0).*((x-c0).^2-Rjj0^2)./((x-c0).^2+Rjj0^2); 
plot(x,V1,'r-',x,norm_nabla_V1,'r--')