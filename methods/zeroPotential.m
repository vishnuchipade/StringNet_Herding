function F=zeroPotential(x)
% ------------------------------
% Author: Vishnu S. Chipade
% Email: vishnuc@umich.edu
% ------------------------------

global rS rO rotSense Q G RO
F(1)=Q*(rS(1)-x(1));
 F(2)= Q*(rS(2)-x(2));
for k=1:size(rO,2)
F(1)=F(1)+(-1)^rotSense(k)*G/(sqrt((rO(1,k)-x(1))^2+(rO(2,k)-x(2))^2)-RO)*(rO(2,k)-x(2))/sqrt((rO(1,k)-x(1))^2+(rO(2,k)-x(2))^2);
 F(2)=F(2)+(-1)^(rotSense(k)+1)*G/(sqrt((rO(1,k)-x(1))^2+(rO(2,k)-x(2))^2)-RO)*(rO(1,k)-x(1))/sqrt((rO(1,k)-x(1))^2+(rO(2,k)-x(2))^2);
end
end