function PlotShortestPath(tanG,Path,tanG_new,figNumber,pathColor,lineType,lineWidth)
% ------------------------------
% Author: Vishnu S. Chipade
% Email: vishnuc@umich.edu
% ------------------------------




global rho_safe
figure(figNumber);
segType=Path.segType;
NS=length(segType);
rV=Path.rV;
rVC=Path.rVC;
%minimum required vertices
hold on
plot(rV(1,1),rV(2,1),'.','color',pathColor)

for i=1:NS
    hold on
    if mod(i,2)==1
        plot(rV(1,i:i+1),rV(2,i:i+1),'color',pathColor,'linestyle',lineType,'linewidth',lineWidth);
    else
        
        ang1=atan2(rV(2,i)-rVC(2,i),rV(1,i)-rVC(1,i));
         if ang1<0
            ang1=ang1+2*pi;
        end
        rotMat=[cos(ang1), sin(ang1);  -sin(ang1), cos(ang1)];
        r2_prime=rotMat*(rV(:,i+1)-rVC(:,i));
        ang2_prime=atan2(r2_prime(2),r2_prime(1));
       
        if segType(i)==1
            if ang2_prime<0
                ang2_prime=ang2_prime+2*pi;
            end
        plot(rVC(1,i)+rho_safe*cos(ang1:(ang2_prime)/50:ang1+ang2_prime),rVC(2,i)+rho_safe*sin(ang1:(ang2_prime)/50:ang1+ang2_prime),'color',pathColor,'linestyle',lineType,'linewidth',lineWidth);
        else
            if ang2_prime>0
                ang2_prime=ang2_prime-2*pi;
            end
         plot(rVC(1,i)+rho_safe*cos(ang1:(ang2_prime)/50:ang1+ang2_prime),rVC(2,i)+rho_safe*sin(ang1:(ang2_prime)/50:ang1+ang2_prime),'color',pathColor,'linestyle',lineType,'linewidth',lineWidth);

        end
            
    end
    hold on
    plot(rV(1,i+1),rV(2,i+1),'.','color',pathColor)
end
end