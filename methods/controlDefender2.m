function [ uD, arr_EDO ] = controlDefender2(XD, XD_des, XD_des_dot, ND )
% ------------------------------
% Author: Vishnu S. Chipade
% Email: vishnuc@umich.edu
% ------------------------------

%%Control action for the defenders to reach their desired positions while avoiding
%obstacles and defenders using vector fields
global rO nO aO bO kr1 kr2 kr0 RD_con u_maxD
global E_m_DO E_bar_DO E_u_DO R_m_DD R_bar_DD R_u_DD
global A_D_O B_D_O C_D_O D_D_O A_D_D B_D_D C_D_D D_D_D
NO=size(rO,2);
for j=1:ND
    F_d=0;  %for sum
    temp2=1; %for product
    temp3=0;
    rD_des=XD_des(1:2,j);
    rD=XD(1:2,j);
    vD=XD(3:4,j);
    vD_des=XD_des(3:4,j);
    norm_drD=norm((rD_des-rD));
    norm_dvD=norm((vD_des-vD));
    % check for nearby obstacles
    if(1)
        for k=1:NO
            n=nO(k);
            a=aO(k);
            b=bO(k);
            EmDO=E_m_DO(j,k);
            EbarDO=E_bar_DO(j,k);
            EuDO=E_u_DO(j,k);
            E_DOk=(abs((rD(1)-rO(1,k))/a))^(2*n)+(abs((rD(2)-rO(2,k))/b))^(2*n)-1;
            E_Ddes_Ok=(abs((rD_des(1)-rO(1,k))/a))^(2*n)+(abs((rD_des(2)-rO(2,k))/b))^(2*n)-1;
            if E_Ddes_Ok<EuDO
                EuDO=E_Ddes_Ok;
                EbarDO=(E_Ddes_Ok-EmDO)/2;
            end
            arr_EDO(j,k)=E_DOk;
            if E_DOk>EmDO &&  E_DOk<EbarDO
                sigma=1;
            elseif  E_DOk>EbarDO &&  E_DOk<EuDO
                sigma=A_D_O(j,k)*E_DOk^3+B_D_O(j,k)*E_DOk^2+C_D_O(j,k)*E_DOk+D_D_O(j,k);
            else
                sigma=0;
            end
            
            %find vector field direction
            betaOS=atan2(rD_des(2)-rO(2,k),rD_des(1)-rO(1,k));  %angle between the desired location and the obstacle
            betaOF=atan2(rD(2)-rO(2,k),rD(1)-rO(1,k));   %angle between the location and the obstacle
            beta_barOS=atan2(b^(2*n)*sign(rD_des(1)-rO(1,k))*(abs(rD_des(1)-rO(1,k)))^(2*n-1),-a^(2*n)*sign(rD_des(2)-rO(2,k))*(abs(rD_des(2)-rO(2,k)))^(2*n-1));  %angle of the tangent at the safe area
          
            beta_bar=atan2(b^(2*n)*sign(rD(1)-rO(1,k))*(abs(rD(1)-rO(1,k)))^(2*n-1),-a^(2*n)*sign(rD(2)-rO(2,k))*(abs(rD(2)-rO(2,k)))^(2*n-1));  %angle of the tangent at the given location
            if beta_bar<0   % get betaOS between [0,2pi]
                beta_bar=beta_bar+2*pi;
            end
            Delta_beta=betaOF-betaOS;
            if Delta_beta<0
                Delta_beta=Delta_beta+2*pi;
            end
            
            Delta_beta_bar=beta_barOS-betaOS;
            if Delta_beta_bar<0
                Delta_beta_bar=Delta_beta_bar+2*pi;
            end
            if Delta_beta<pi
                beta_bar1=beta_bar-Delta_beta_bar+Delta_beta/pi*(Delta_beta_bar-pi);
            else
                beta_bar1=beta_bar-Delta_beta_bar/pi*(Delta_beta-pi);
            end
                       
            F_d = F_d + sigma*[cos(beta_bar1);sin(beta_bar1)];
            temp2=temp2*(1-sigma);
            temp3=temp3+sigma;
        end
    end
    
    %check for nearby defenders
    if(1)
        for l=1:ND
            if l~=j
                R_ADjDl=norm(rD-XD(1:2,l));
                if R_ADjDl>R_m_DD && R_ADjDl<R_bar_DD
                    sigma=1;
                elseif R_ADjDl>R_bar_DD && R_ADjDl<R_u_DD
                    sigma=A_D_D*R_ADjDl^3+B_D_D*R_ADjDl^2+C_D_D*R_ADjDl+D_D_D;
                else
                    sigma=0;
                end
                F_d = F_d+ sigma*(rD-XD(1:2,l))/R_ADjDl;
                temp2=temp2*(1-sigma);
                temp3=temp3+sigma;
            end
        end
        %temp3
    end
    %for attraction toward the desired locations
    if norm_drD~=0
        F_d = F_d + temp2*(rD_des-rD)/norm_drD;%/Rc^2;
    end
    if norm(F_d)~=0
        if norm_drD>RD_con
            vD_des0=kr0*tanh(norm_drD)*F_d/norm(F_d);
            %norm_dvD=norm((vD_des0-vD));
            uD_des0=kr0*tanh(norm_dvD)*F_d/norm(F_d)+vD_des0;
        else
            vD_des0=kr1*norm_drD^kr2*F_d/norm(F_d);
            % norm_dvD=norm((vD_des0-vD));
            uD_des0=0.5*kr1*norm_dvD^kr2*F_d/norm(F_d)+vD_des0;
        end
        %F_d
    else
        uD_des0=1;
    end
%    vD_des0=kr0*tanh(norm_drD)*F_d/norm(F_d);
%             uD_des0=kr0*tanh(norm_dvD)*F_d/norm(F_d)+vD_des0;
    
    if temp3==0
        uD(:,j)=XD_des_dot(3:4,j)+uD_des0;        
    else
       uD(:,j)=uD_des0;
        %2
    end
    norm_uD=norm(uD);
    if norm_uD>u_maxD(j)
       uD=uD*u_maxD(j)/norm_uD;
    end
    
end
end
