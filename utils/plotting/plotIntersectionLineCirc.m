function plotIntersectionLineCirc(r1,r2,r21,r22,rVC2,segType2,interSec)
% ------------------------------
% Author: Vishnu S. Chipade
% Email: vishnuc@umich.edu
% ------------------------------

global rho_safe rho_D


circCos=cos(0:pi/100:2*pi);
circSin=sin(0:pi/100:2*pi);
fontSize=18;

figure
hold all
axis equal
%Plot circle in light color
plot(rVC2(1)+rho_safe*circCos,rVC2(2)+rho_safe*circSin,'-.','color',[0,0,0.5])
%Plot circular segment
ang1=atan2(r21(2)-rVC2(2),r21(1)-rVC2(1));
if ang1<0
    ang1=ang1+2*pi;
end
rotMat=[cos(ang1), sin(ang1);  -sin(ang1), cos(ang1)];
r2_prime=rotMat*(r22(:)-rVC2(:));
ang2_prime=atan2(r2_prime(2),r2_prime(1));
pathColor=[0,0,1];
if segType2==1
    if ang2_prime<0
        ang2_prime=ang2_prime+2*pi;
    end
    plot(rVC2(1)+rho_safe*cos(ang1:(ang2_prime)/50:ang1+ang2_prime),rVC2(2)+rho_safe*sin(ang1:(ang2_prime)/50:ang1+ang2_prime),'color',pathColor);
    plot(rVC2(1)+(rho_safe-rho_D)*cos(ang1:(ang2_prime)/50:ang1+ang2_prime),rVC2(2)+(rho_safe-rho_D)*sin(ang1:(ang2_prime)/50:ang1+ang2_prime),'--','color',pathColor);
    plot(rVC2(1)+(rho_safe+rho_D)*cos(ang1:(ang2_prime)/50:ang1+ang2_prime),rVC2(2)+(rho_safe+rho_D)*sin(ang1:(ang2_prime)/50:ang1+ang2_prime),'--','color',pathColor);

else
    if ang2_prime>0
        ang2_prime=ang2_prime-2*pi;
    end
    plot(rVC2(1)+rho_safe*cos(ang1:(ang2_prime)/50:ang1+ang2_prime),rVC2(2)+rho_safe*sin(ang1:(ang2_prime)/50:ang1+ang2_prime),'color',pathColor);
    plot(rVC2(1)+(rho_safe-rho_D)*cos(ang1:(ang2_prime)/50:ang1+ang2_prime),rVC2(2)+(rho_safe-rho_D)*sin(ang1:(ang2_prime)/50:ang1+ang2_prime),'--','color',pathColor);
    plot(rVC2(1)+(rho_safe+rho_D)*cos(ang1:(ang2_prime)/50:ang1+ang2_prime),rVC2(2)+(rho_safe+rho_D)*sin(ang1:(ang2_prime)/50:ang1+ang2_prime),'--','color',pathColor);

end
text(r21(1),r21(2),'$r_{21}$','fontsize',fontSize)
text(r22(1),r22(2),'$r_{22}$','fontsize',fontSize)

%Plot straight line
plot([r1(1),r2(1)],[r1(2),r2(2)],'r')
text(r1(1),r1(2),'$r_1$','fontsize',fontSize)
text(r2(1),r2(2),'$r_2$','fontsize',fontSize)



%plot intersection points on the line
p1=interSec.Pos1(:,1);
p2=interSec.Pos1(:,2);
plot(p1(1)+rho_D*circCos,p1(2)+rho_D*circSin,'r--')
plot(p2(1)+rho_D*circCos,p2(2)+rho_D*circSin,'r--')
plot([p1(1),p2(1)],[p1(2),p2(2)],'--');

%plot intersection points on the arc
p1=interSec.Pos2(:,1);
p2=interSec.Pos2(:,2);
plot(p1(1)+rho_D*circCos,p1(2)+rho_D*circSin,'b--')
plot(p2(1)+rho_D*circCos,p2(2)+rho_D*circSin,'b--')
end