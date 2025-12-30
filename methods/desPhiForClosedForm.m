function phi=desPhiForClosedForm(rAc)
% ------------------------------
% Author: Vishnu S. Chipade
% Email: vishnuc@umich.edu
% ------------------------------

global rP obs A_Ac_Oc B_Ac_Oc C_Ac_Oc D_Ac_Oc R_bar_AcOc R_u_AcOc
NO=obs.NO;
rCO2=obs.rCO2;
sigmaSum=0;
sigmaProd=1; %for product

for k=1:NO
    RAcOc=norm(rAc-rCO2(:,k));
    if RAcOc<R_bar_AcOc(k)
        sigma=1;
    elseif RAcOc>=R_bar_AcOc(k) && RAcOc<=R_u_AcOc(k)
        sigma=A_Ac_Oc(k)*RAcOc^3+B_Ac_Oc(k)*RAcOc^2+C_Ac_Oc(k)*RAcOc+D_Ac_Oc(k);
    else
        sigma=0;
    end
phiAcOc=atan2(rAc(2)-rCO2(2,k),rAc(1)-rCO2(1,k));
if phiAcOc<0
phiAcOc=phiAcOc+2*pi;
end
    sigmaSum=sigmaSum+sigma*phiAcOc;
    sigmaProd=sigmaProd*(1-sigma);    
end
phiAcP=atan2(rAc(2)-rP(2),rAc(1)-rP(1));
if phiAcP<0
phiAcP=phiAcP+2*pi;
end

phi=sigmaProd*phiAcP+sigmaSum;
end