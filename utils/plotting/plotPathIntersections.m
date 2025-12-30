function [interSec]=plotPathIntersections(Path1,Path2,tanG,tanG_new1,tanG_new2,figNumber,pathColors,lineWidth)
% ------------------------------
% Author: Vishnu S. Chipade
% Email: vishnuc@umich.edu
% ------------------------------


global rho_D

colors=distinguishable_colors(20);

plotRec=1;  %if you want to plot rectangles to visualize

circCos=cos(0:pi/50:2*pi);
circSin=sin(0:pi/50:2*pi);

[interSec]=pathIntersections(Path1,Path2);

if sum(sum(interSec.Flag))>0
    %Plot the plots
    %pathColor='r';
    fontSize=18;
    figure(figNumber)
    hold on;
    
    PlotShortestPath(tanG,Path1,tanG_new1,figNumber,pathColors{1},'-',lineWidth)
    PlotShortestPath(tanG,Path2,tanG_new2,figNumber,pathColors{2},'-',lineWidth)
    for i=1:size(interSec.rV11,2)
        %On path 1
        rP1=interSec.rV11(:,i);
        rP2=interSec.rV12(:,i);
        plot(rP1(1)+rho_D*circCos,rP1(2)+rho_D*circSin,'--','color',pathColors{1}(1:3))
        plot(rP2(1)+rho_D*circCos,rP2(2)+rho_D*circSin,'--','color',pathColors{1}(1:3))
        if plotRec==1
        rPmid=(rP1+rP2)/2;
        Pbar=norm(rP2-rP1);
        theta1=atan2(rP2(2)-rP1(2),rP2(1)-rP1(1));
        rec1=[rPmid(1)*ones(1,5);rPmid(2)*ones(1,5)]+[cos(theta1),-sin(theta1);sin(theta1),cos(theta1)]*([-Pbar/2,Pbar/2,Pbar/2,-Pbar/2,-Pbar/2;-rho_D,-rho_D,rho_D,rho_D,-rho_D])              
        plot(rec1(1,:),rec1(2,:),'-','color',pathColors{1}(1:3))  %swapping rectangle
        end
        text(rP1(1),rP1(2),num2str(i),'fontsize',fontSize,'color',colors(i,:))
    end
    for j=1:size(interSec.rV21,2)
        %On path 2
        rP1=interSec.rV21(:,j);
        rP2=interSec.rV22(:,j);        
        plot(rP1(1)+rho_D*circCos,rP1(2)+rho_D*circSin,'color',pathColors{2}(1:3))
        plot(rP2(1)+rho_D*circCos,rP2(2)+rho_D*circSin,'color',pathColors{2}(1:3))
        if plotRec==1
        rPmid=(rP1+rP2)/2;
        Pbar=norm(rP2-rP1);
        theta1=atan2(rP2(2)-rP1(2),rP2(1)-rP1(1));
        rec1=[rPmid(1)*ones(1,5);rPmid(2)*ones(1,5)]+[cos(theta1),-sin(theta1);sin(theta1),cos(theta1)]*([-Pbar/2,Pbar/2,Pbar/2,-Pbar/2,-Pbar/2;-rho_D,-rho_D,rho_D,rho_D,-rho_D])
        plot(rec1(1,:),rec1(2,:),'-','color',pathColors{2}(1:3))  %swapping rectangle
        end
        %text(rP1(1),rP1(2),num2str(j),'fontsize',fontSize,'color',colors(j,:))
    end
   
end