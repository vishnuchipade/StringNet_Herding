
%This script is used to plot the maximum and minimum values of the angle

% ------------------------------
% Author: Vishnu S. Chipade
% Email: vishnuc@umich.edu
% ------------------------------

%between two vector fields (attractive and the rectangular repulsive) in the blending region
if(0)
    figure
    sig=0:0.01:1;
    plot(sig,sig.*(1-sig))
end

plot_delBeta=0;
for bt=1:180
    Bt(bt)=bt/2;
    beta_rs=Bt(bt)*pi/180;
    rS=5*[cos(beta_rs);sin(beta_rs)];
    n=1.2;
    a=2;
    b=1.5;
    E=(rS(1)/a)^(2*n)+(rS(2)/b)^(2*n)-1;
    const1=a*(1+E)^(0.5/n);
    const2=b*(1+E)^(0.5/n);
    beta_rs=atan2(rS(2),rS(1));
    t=beta_rs-pi:pi/100:beta_rs-pi/100;
    %t=0:pi/100:pi/2;
    beta_bar_rs=atan2(cos(beta_rs)^2,-sin(beta_rs)^2)
    for i=1:length(t)
        beta_bar=atan2(sign(cos(t(i)))*cos(t(i))^2,-sign(sin(t(i)))*sin(t(i))^2)
        %     if beta_bar<0
        %         beta_bar=beta_bar+2*pi;
        %     end
        dbeta(i)=beta_bar-t(i);
        if t(i)~=beta_rs
            delBeta(i)=atan2(rS(2)-const2*sign(sin(t(i)))*abs(sin(t(i)))^(1/n),rS(1)-const1*sign(cos(t(i)))*abs(cos(t(i)))^(1/n))-beta_bar+(beta_bar_rs-beta_rs)/pi*(pi+t(i)-beta_rs);
        else
            delBeta(i)=(beta_bar_rs-beta_rs)/pi*(pi+t(i)-beta_rs);
        end
        x(i)=const1*sign(cos(t(i)))*abs(cos(t(i)))^(1/n);
        y(i)=const2*sign(sin(t(i)))*abs(sin(t(i)))^(1/n);
    end
    min_delBeta(bt)=min(delBeta);
    max_delBeta(bt)=max(delBeta);
    if(plot_delBeta==1)
    figure
    plot(t,dbeta)
    hold on;
    plot(t,delBeta,'r')
    xlabel('$\beta$')
    end
    
end
%Plot the minimum and the maximum values
figure
plot(Bt,min_delBeta,Bt,max_delBeta,'r')
ylabel('$\partial{\bar \beta_{ok}}$')
xlabel('$\beta_{ok}(\textbf{r}_s)$')
legend('$\min(\partial{\bar \beta_{ok}})$','$\max(\partial{\bar \beta_{ok}})$')

if(0)
    f0=@(x,y) abs((x/a))^(2*n)+abs((y/b))^(2*n)-1-E
    figure
    ezplot(f0,[-10,10])
    hold on
    plot(x,y,'r')
end

f=@(x) 2*cos(x)*sin(x)-cos(x)^4-sin(x)^4
beta0=fsolve(f,0.9)
dbeta0=atan2(cos(beta0)^2,-sin(beta0)^2)-beta0;

a=beta_rs; b=(beta_bar_rs-beta_rs)/pi;
f1=@(x) -(sec(x)^2)/((tan(a) - tan(x))^2 + 1) - b - (4*sin(2*x))/(cos(4*x) + 3)
beta0=fsolve(f1,beta_rs-pi)
delBeta0=atan2(tan(beta_rs)-tan(beta0),1)-atan2(cos(beta0)^2,-sin(beta0)^2)+(beta_bar_rs-beta_rs)/pi*(pi+beta_rs-beta0)
