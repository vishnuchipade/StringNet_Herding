% ------------------------------
% Author: Vishnu S. Chipade
% Email: vishnuc@umich.edu
% ------------------------------

syms Del p p_dot p_ddot pp pp_dot  g g_dot g_ddot FAO FAO_dot FAO_ddot

Delta_barS=1;
F_bar_AO=0.9;
f=@(p,p1,g,Del,FAO) Del*(sin(pp)-tan(p)*cos(pp))-FAO*(tan(p)*cos(g)-sin(g));

dfdp=diff(f,p);
dfdpp=diff(f,pp);
dfdg=diff(f,g);
dfdFAO=diff(f,FAO);
d2fdp2=diff(f,p,2);
d2fdpp2=diff(f,pp,2);
d2fdg2=diff(f,g,2);
d2fdFAO2=diff(f,FAO,2);
psi_prime_dot0=-(dfdp*p_dot+dfdg*g_dot+dfdFAO*FAO_dot)/dfdpp;
psi_prime_dot_fun=matlabFunction(psi_prime_dot0,'Vars',[Del,FAO, FAO_dot, pp,p,p_dot,g,g_dot]);
psi_prime_ddot0=-(d2fdp2*p_dot^2+dfdp*p_ddot+d2fdg2*g_dot^2+dfdg*g_ddot+d2fdpp2*pp_dot^2+d2fdFAO2*FAO_dot^2+dfdFAO*FAO_ddot)/dfdpp;
psi_prime_ddot_fun=matlabFunction(psi_prime_ddot0,'Vars',[Del,FAO, FAO_dot, FAO_ddot, pp,pp_dot,p,p_dot,p_ddot,g,g_dot,g_ddot]);


%For sigma product for the formation 
global sigmaProd_dot_fun sigmaProd_ddot_fun
sig = sym('sig',[1,NO]);
sig_dot=sym('sig_dot',[1,NO]);
sig_ddot=sym('sig_ddot',[1,NO]);
sigProd=1;
for k=1:NO
    sigProd=sigProd*(1-sig(k));
end

for k=1:NO
    dsigProd(k)=diff(sigProd,sig(k));
    d2sigProd(k)=diff(sigProd,sig(k),2);
end
if exist('dsigProd','var')
sigProd_dot=sum(dsigProd.*sig_dot);
sigmaProd_dot_fun=matlabFunction(sigProd_dot,'Vars',[sig sig_dot]);

sigProd_ddot=sum(d2sigProd.*sig_dot.^2+dsigProd.*sig_ddot);
sigmaProd_ddot_fun=matlabFunction(sigProd_ddot,'Vars',[sig sig_dot sig_ddot]);
end

%Sigma product for the attacker
global sigmaProdA_dot_fun 
NODA=NO+ND+NA;
sigA = sym('sigA',[1,NODA]);
sigA_dot=sym('sigA_dot',[1,NODA]);
sigProdA=1;
for k=1:NODA
    sigProdA=sigProdA*(1-sigA(k));
end

for k=1:NODA
     dsigProdA(k)=diff(sigProdA,sigA(k));
%     d2sigProd(k)=diff(sigProd,sig(k),2);
end
sigProdA_dot=sum(dsigProdA.*sigA_dot);
sigmaProdA_dot_fun=matlabFunction(sigProdA_dot,'Vars',[sigA sigA_dot]);

%Sigma product for the attacker
global sigmaProdD_dot_fun 
NOD=NO+ND-1;
sigD = sym('sigD',[1,NOD]);
sigD_dot=sym('sigD_dot',[1,NOD]);
sigProdD=1;
for k=1:NOD
    sigProdD=sigProdD*(1-sigD(k));
end

for k=1:NOD
     dsigProdD(k)=diff(sigProdD,sigD(k));
%     d2sigProd(k)=diff(sigProd,sig(k),2);
end
sigProdD_dot=sum(dsigProdD.*sigD_dot);
sigmaProdD_dot_fun=matlabFunction(sigProdD_dot,'Vars',[sigD sigD_dot]);

global d2Beta_barF_dBeta2_fun
syms beta n0 b0 a0 d2Beta_barF_dBeta2_fun
beta_barF=atan2(b0^(2*n0)*cos(beta)^(2*n0-1),-a0^(2*n0)*sin(beta)^(2*n0-1));
dBeta_barF_dBeta=diff(beta_barF, beta);
d2Beta_barF_dBeta2=diff(dBeta_barF_dBeta, beta);
d2Beta_barF_dBeta2_fun=matlabFunction(d2Beta_barF_dBeta2,'Vars',[beta n0 a0 b0]);
% sigProd_ddot=sum(d2sigProd.*sig_dot.^2+dsigProd.*sig_ddot);
% sigmaProd_ddot_fun=matlabFunction(sigProd_ddot,'Vars',[sig sig_dot sig_ddot]);
