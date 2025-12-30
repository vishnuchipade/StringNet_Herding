% ------------------------------
% Author: Vishnu S. Chipade
% Email: vishnuc@umich.edu
% ------------------------------

a=2;
b=1;
R=2;
x1=a*cos(0:pi/100:2*pi);
y1=b*sin(0:pi/100:2*pi);
x2=(a+R)*cos(0:pi/100:2*pi);
y2=(b+R)*sin(0:pi/100:2*pi);

for i=1:length(x1);
    if y1(i)>=0
        m1=-2*b*x1(i)/sqrt(a^2-x1(i)^2);
        m2=-1/m1;
        theta1=atan(m2);
        if theta1<0
            theta1=theta1+pi;
        end
        x3(i)=x1(i)+R*cos(theta1);
        y3(i)=y1(i)+R*sin(theta1);
    else
         m1=2*b*x1(i)/sqrt(a^2-x1(i)^2);
        m2=-1/m1;
        theta1=atan(m2);
        if theta1>0
            theta1=theta1+pi;
        end
        x3(i)=x1(i)+R*cos(theta1);
        y3(i)=y1(i)+R*sin(theta1);
    end
end
    
figure
plot(x1,y1,x2,y2,x3,y3)
