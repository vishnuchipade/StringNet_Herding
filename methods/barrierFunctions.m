%Plots barrier functions 

% ------------------------------
% Author: Vishnu S. Chipade
% Email: vishnuc@umich.edu
% ------------------------------


rP=[0,0]';
rO=[-10,0;0,-10]'%;0,10;10,0;30,0;40,0]';
NO=size(rO,2);
kappa_AO=60;
eps_AO=0.1;
syms x y
f1 = x-rP(1);
f2 = y-rP(2);
for k=1:NO
f1= f1-2*(x-rO(1,k))*kappa_AO^2/((x-rO(1,k))^2+(y-rO(2,k))^2+eps_AO^2)^2;%-2*(x-rO(1,2))*kappa_AO^2/((x-rO(1,2))^2+(y-rO(2,2))^2+eps_AO^2)^2 ...
    %-2*(x-rO(1,3))*kappa_AO^2/((x-rO(1,3))^2+(y-rO(2,3))^2+eps_AO^2)^2;
f2= f2-2*(y-rO(2,k))*kappa_AO^2/((x-rO(1,k))^2+(y-rO(2,k))^2+eps_AO^2)^2;%-2*(y-rO(2,2))*kappa_AO^2/((x-rO(1,2))^2+(y-rO(2,2))^2+eps_AO^2)^2 ...
    %-2*(y-rO(2,3))*kappa_AO^2/((x-rO(1,3))^2+(y-rO(2,3))^2+eps_AO^2)^2;
end
 f=@(x,y) x-rP(1)-2*(x-rO(1,1))*kappa_AO^2/((x-rO(1,1))^2+(y-rO(2,1))^2+eps_AO^2)^2-2*(x-rO(1,2))*kappa_AO^2/((x-rO(1,2))^2+(y-rO(2,2))^2+eps_AO^2)^2 ...
     == y-rP(2)-2*(y-rO(2,1))*kappa_AO^2/((x-rO(1,1))^2+(y-rO(2,1))^2+eps_AO^2)^2-2*(y-rO(2,2))*kappa_AO^2/((x-rO(1,2))^2+(y-rO(2,2))^2+eps_AO^2)^2;

%  f1=@(x,y) f1;
%  f2=@(x,y) f2;
 % f=solve({f1,f2})
%ezplot(f)
%f1=matlabfunction(f1,[x,y])
% f1=matlabFunction(f1,'Vars',[x,y]);
% f2=matlabFunction(f2,'Vars',[x,y]);
figure
plot(rP(1),rP(2),'b*')
hold on
for k=1:NO
plot(rO(1,k),rO(2,k),'ko')
hold on
end
xl=50;
ez1=ezplot(f1,[-xl,xl]);
set(ez1,'color',[1 0 0])
hold on
ezplot(f2,[-xl,xl])

if(1)
[X,Y]=meshgrid(-xl:.1:xl,-xl:.1:xl);
 clear x y
 f10=@(x,y) f1;
 f20=@(x,y) f2;
for i=1:size(X,1)
    for j=1:size(X,2)
        x=X(i,j);  y=Y(i,j);
        F1(i,j)=-f10(X(i,j),Y(i,j));
        F2(i,j)=-f20(X(i,j),Y(i,j));
    end
end
figure
quiver(X,Y,F1,F2,5)
hold on
contour(X,Y,F1+F2,0)
end
