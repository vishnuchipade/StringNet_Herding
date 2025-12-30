function interSec=interSecCircCirc(r11,r12,rotMat,ang12_prime,segType1,r21,r22,segType2,rVC)
% ------------------------------
% Author: Vishnu S. Chipade
% Email: vishnuc@umich.edu
% ------------------------------

%This function assumes that the segments lie on the same circular arc and
%finds the common portion on these segments where the robots can collide
global dthetai rho_D

r21_prime=rotMat*(r21-rVC);   %Rotate relative to the center of arc
ang21_prime=atan2(r21_prime(2),r21_prime(1));
r22_prime=rotMat*(r22-rVC);
ang22_prime=atan2(r22_prime(2),r22_prime(1));

interSec.angfi=[];
angfi11=[];
if segType1==1  %if the first arc is anticlockwise
    angcheck1= 0<=ang21_prime && ang21_prime<=ang12_prime;
    angcheck2= 0<=ang22_prime && ang22_prime<=ang12_prime;
    if segType2==1  %if the second arc is anticlockwise
        if angcheck1 && angcheck2
            angfi21=ang21_prime;
            angfi22=ang22_prime;
            angfi11=max(0,ang21_prime-dthetai);
            angfi12=min(ang12_prime,ang22_prime+dthetai);
        elseif ang22_prime>ang12_prime && angcheck1
            angfi21=ang21_prime;
            angfi22=min(ang22_prime,ang12_prime+dthetai);
            angfi11=max(0,ang21_prime-dthetai);
            angfi12=ang12_prime;
        elseif ang21_prime<0 && angcheck2
            angfi21=max(ang21_prime,-dthetai);
            angfi22=ang22_prime;
            angfi11=0; %ang11_prime=0
            angfi12=min(ang12_prime,ang22_prime+dthetai);
        elseif ang21_prime<0 && ang22_prime>ang12_prime
            angfi21=max(ang21_prime,-dthetai);
            angfi22=min(ang22_prime,ang12_prime+dthetai);
            angfi11=0;
            angfi12=ang12_prime;
        elseif ang21_prime<0 && ang22_prime<0
            if norm(r11-r22)<2*rho_D
                angfi21=max(ang21_prime,ang22_prime-dthetai);
                angfi22=ang22_prime;
                angfi11=0;
                angfi12=min(ang12_prime,0+dthetai);
            end
        elseif ang21_prime>ang12_prime && ang22_prime>ang12_prime
            if norm(r12-r21)<2*rho_D
                angfi21=ang21_prime;
                angfi22=min(ang22_prime,ang21_prime+dthetai);
                angfi11=max(0,ang12_prime-dthetai);
                angfi12=ang12_prime;
            end
        end
    elseif segType2==2   %if the second arc is clockwise
        if angcheck1 && angcheck2
            angfi21=ang21_prime;
            angfi22=ang22_prime;
            angfi11=max(0,ang22_prime-dthetai);
            angfi12=min(ang12_prime,ang21_prime+dthetai);
        elseif ang21_prime>ang12_prime && angcheck2
            angfi22=ang22_prime;
            angfi21=min(ang21_prime,ang12_prime+dthetai);
            angfi11=max(0,ang22_prime-dthetai);
            angfi12=ang12_prime;
        elseif ang22_prime<0 && angcheck1
            angfi22=max(ang22_prime,-dthetai);
            angfi21=ang21_prime;
            angfi11=0; %ang11_prime=0
            angfi12=min(ang12_prime,ang21_prime+dthetai);
        elseif ang22_prime<0 && ang21_prime>ang12_prime
            angfi22=max(ang22_prime,-dthetai);
            angfi21=min(ang21_prime,ang12_prime+dthetai);
            angfi11=0;
            angfi12=ang12_prime;
        elseif ang21_prime<0 && ang22_prime<0
            if norm(r11-r21)<2*rho_D
                angfi22=max(ang22_prime,ang21_prime-dthetai);
                angfi21=ang21_prime;
                angfi11=0;
                angfi12=min(ang12_prime,0+dthetai);
            end
        elseif ang21_prime>ang12_prime && ang22_prime>ang12_prime
            if norm(r12-r22)<2*rho_D
                angfi22=ang22_prime;
                angfi21=min(ang21_prime,ang22_prime+dthetai);
                angfi11=max(0,ang12_prime-dthetai);
                angfi12=ang12_prime;
            end
        end
    end
elseif segType1==2   %if the first arc is clockwise
    angcheck1= 0>=ang21_prime && ang21_prime>=ang12_prime;
    angcheck2= 0>=ang22_prime && ang22_prime>=ang12_prime;
    if  segType2==2   %if the second arc is clockwise
        if angcheck1 && angcheck2
            angfi21=ang21_prime;
            angfi22=ang22_prime;
            angfi11=min(0,ang21_prime+dthetai);
            angfi12=max(ang12_prime,ang22_prime-dthetai);
        elseif ang22_prime<ang12_prime && angcheck1
            angfi21=ang21_prime;
            angfi22=max(ang22_prime,ang12_prime-dthetai);
            angfi11=min(0,ang21_prime+dthetai);
            angfi12=ang12_prime;
        elseif ang21_prime>0 && angcheck2
            angfi21=min(ang21_prime,dthetai);
            angfi22=ang22_prime;
            angfi11=0; %ang11_prime=0
            angfi12=max(ang12_prime,ang22_prime-dthetai);
        elseif ang21_prime>0 && ang22_prime<ang12_prime
            angfi21=min(ang21_prime,dthetai);
            angfi22=max(ang22_prime,ang12_prime-dthetai);
            angfi11=0;
            angfi12=ang12_prime;
        elseif ang21_prime>0 && ang22_prime>0
            if norm(r11-r22)<2*rho_D
                angfi21=min(ang21_prime,ang22_prime+dthetai);
                angfi22=ang22_prime;
                angfi11=0;
                angfi12=max(ang12_prime,0-dthetai);
            end
        elseif ang21_prime<ang12_prime && ang22_prime<ang12_prime
            if norm(r12-r21)<2*rho_D
                angfi21=ang21_prime;
                angfi22=max(ang22_prime,ang21_prime-dthetai);
                angfi11=min(0,ang12_prime+dthetai);
                angfi12=ang12_prime;
            end
        end
    elseif segType2==1   %if the second arc is anticlockwise        
        if angcheck1 && angcheck2
            angfi21=ang21_prime;
            angfi22=ang22_prime;
            angfi11=min(0,ang22_prime+dthetai);
            angfi12=max(ang12_prime,ang21_prime-dthetai);
        elseif ang21_prime<ang12_prime && angcheck2
            angfi22=ang22_prime;
            angfi21=max(ang21_prime,ang12_prime-dthetai);
            angfi11=min(0,ang22_prime+dthetai);
            angfi12=ang12_prime;
        elseif ang22_prime>0 && angcheck1
            angfi22=min(ang22_prime,dthetai);
            angfi21=ang21_prime;
            angfi11=0; %ang11_prime=0
            angfi12=max(ang12_prime,ang21_prime-dthetai);
        elseif ang22_prime>0 && ang21_prime<ang12_prime
            angfi22=min(ang22_prime,dthetai);
            angfi21=max(ang21_prime,ang12_prime-dthetai);
            angfi11=0;
            angfi12=ang12_prime;
        elseif ang21_prime>0 && ang22_prime>0
            if norm(r11-r21)<2*rho_D
                angfi22=min(ang22_prime,ang21_prime+dthetai);
                angfi21=ang21_prime;
                angfi11=0;
                angfi12=max(ang12_prime,0-dthetai);
            end
        elseif ang21_prime<ang12_prime && ang22_prime<ang12_prime
            if norm(r12-r22)<2*rho_D
                angfi22=ang22_prime;
                angfi21=max(ang21_prime,ang22_prime-dthetai);
                angfi11=min(0,ang12_prime+dthetai);
                angfi12=ang12_prime;
            end
        end 
    end
end
if ~isempty(angfi11)
    interSec.angfi=[angfi11,angfi12,angfi21,angfi22];
end
end
