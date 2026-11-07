function dxdt = modelLIH_RA2019(X,p)

%Unpack sate vector and convert from uM -> M
H      = X(1)*6.02e17;  %umol -> N
pH_L     = X(2);  %unitless
Ca_T   = X(3)/1e6; %uM -> M
Ca_F = p.r*Ca_T;  %M

Ca_C   = X(4)/1e6;  %uM -> M
Ca_md = X(5)/1e6;   %uM -> M
NAADP = X(6);        %uM

delta_pH    = p.pH_C-pH_L;  %unitless


%V-ATPase performance
V_ATPASE    = lys.find_VATP_rate(p.v_flux,-40,pH_L);  %N/ms/pump
J_VATPASE   = p.N_VATP*real(V_ATPASE)*0.015;  %N/ms


%CAX Antiporter {H out, Ca in}
CAX_mu   = p.psi_cax + p.RTF*(p.CAX_H*2.3*delta_pH + p.CAX_Ca/2*log(Ca_F/Ca_C));  %N/s/antiporter
J_CAX    = p.N_CAX*CAX_mu/1000; %N/ms
J_CAX_Ca = (p.CAX_Ca*J_CAX / (p.NA*p.init_V)) * 1e6; %uM/ms for lysosome


%Passive flux
J_H = p.v_H * (10^-(p.pH_C) - 10^(-pH_L)) * p.init_V * p.NA; %N/ms
J_Ca = p.v_Ca * (Ca_C - Ca_F) * 1e6;  %uM/ms


%TPC Channel - Generic model, not based on experimental data
f_NAADP = (NAADP^p.n_TPC/(p.Ka_TPC^p.n_TPC + NAADP^p.n_TPC)); %unitless
J_TPC = p.OCaR*p.v_TPC*f_NAADP*(Ca_md-Ca_F) * 1e6;  %uM/ms


%Time Dependent Quantities
dNHdt   = J_H + (J_VATPASE) - (J_CAX*p.CAX_H);  %N
dHdt  = (dNHdt/(p.NA*p.init_V))*1e6; %uM
dNHdt_umol = dNHdt/6.02e17; %umol

dpHdt   = (-dNHdt/(p.init_V*p.NA))/p.beta_pH; %unitless

dCaTdt = J_Ca + (J_CAX_Ca) + J_TPC;   %uM


J_uM = [J_Ca; J_CAX_Ca; J_TPC; (J_VATPASE/(p.NA*p.init_V))*1e6];
%J_uM = [0; 0; 0; 0];

%OUTPUT
dxdt = [dNHdt_umol; dpHdt; dCaTdt; J_uM];
%dxdt = [0;0;0;J_uM];
