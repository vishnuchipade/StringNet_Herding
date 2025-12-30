function angc=angInterSecCircCirc(rc1,rhoc1,rc2,rhoc2,rotMat)
% ------------------------------
% Author: Vishnu S. Chipade
% Email: vishnuc@umich.edu
% ------------------------------

%this function calcualtes intersection points between two circles and
%their angles measured with respect to line passing through a line joining
%center rc1 and first point on its arc (rotMat gives the corresponding
%rotation matrix for calculations

%intersection points between two circles
d=norm(rc2-rc1);  %Distance between two centers (P0=rc2,P1=rc1)
a=(rhoc2^2-(rhoc1)^2+d^2)/(2*d);
h=sqrt(rhoc2^2-a^2);
P2=rc2(:)+a*(rc1-rc2(:))/d;
if d<rhoc1+rhoc2
%first intersection point
x_int=P2(1)+h*(rc1(2)-rc2(2))/d;
y_int=P2(2)-h*(rc1(1)-rc2(1))/d;
ri0=rotMat*([x_int,y_int]'-rc2(:));
angc(1)=atan2(ri0(2),ri0(1));

%second intersection point
x_int=P2(1)-h*(rc1(2)-rc2(2))/d;
y_int=P2(2)+h*(rc1(1)-rc2(1))/d;
ri0=rotMat*([x_int,y_int]'-rc2(:));
angc(2)=atan2(ri0(2),ri0(1));
else
    angc=[];
end


end