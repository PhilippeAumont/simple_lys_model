function generate_parameters()
%This file creates the file containing all parameters.

N_Lys = 300;            %Number of lysosomes per cell

p.A_cap = 1.534e-4;   %cm2
p.V_cell = 38.00e-6;  %uL
p.V_cyt = 25.84e-6;   %uL
p.V_JSR = 0.12e-6;    %uL
p.V_NSR = 2.098e-6;   %uL
p.V_ss = 1.485e-9;    %uL
p.V_cav = 7.600e-07;  %uL (2% of cell volume)
p.V_ecav = 1.520e-06; %uL (4% of cell volume)

p.V_md = 1.0367e-13*N_Lys; %uL (cylinder with r=100nm and h=3.3nm)
p.init_V = 1.64636595e-16*N_Lys;    %L
p.V_lys = p.init_V * 1e6; %uL

% %Extracellular Ion concentrations
p.K_o = 5400;      %uM
p.Na_o = 140000;   %uM
p.Ca_o = 1800;     %uM

% %SR Parameters
p.v1 = 4.5;          %1/ms %From RyR Module Parameters
p.v2 = 1.74e-5;    %1/ms
p.v3 = 0.306;       %1/ms %From PLB module
p.t_tr = 20;         %ms
p.t_xfer = 8;        %ms

% %Calmodulin and Calsequestrin
p.CMDN_tot = 50.0;      %uM
p.CSQN_tot = 15000.0;  %uM
p.Km_CMDN = 0.238;      %uM
p.Km_CSQN = 800.0;      %uM

% %Membrane current parameters
p.C_m = 1.0;           %uF/cm2
p.F = 96.5;            %C/mmol
p.T = 298;             %K
p.R = 8.314;           %J/mol*K
p.k_NaCa = 275;        %pA/pF
p.Km_Na = 87500;      %uM
p.Km_Ca = 1380;        %uM
p.k_sat = 0.27;
p.n = 0.35;
p.I_max_pCa = 0.051;
p.Km_pCa = 0.5;        %uM
p.G_Cab = 0.000284;   %mS/uF
p.G_Nab = 0.0063;     %mS/uF
p.G_Kss = 0.0611;     %mS/uF
p.G_Ks = 0.00575;     %mS/uF
p.G_Kr = 0.078;        %mS/uF
p.k_f = 0.023761;     %1/ms
p.k_b = 0.036778;     %1/ms
p.G_ClCa = 10.0;       %mS/uF
p.Km_Cl = 10.0;        %uM
p.E_Cl = -40;           %mV

%%Lys parameters
p.CAX_Ca = 1;
p.CAX_H = 3;
p.F_Lys = 96485;      %C/mol
p.NA = 6.02e23;       %u/mol
p.P_Ca = 1.49e-10;   %cm/ms
p.P_Cl = 1.2e-8;      %cm/ms
p.P_H = 6e-8;          %cm/ms
p.P_K = 7.1e-10;      %cm/ms
p.P_Na = 9.6e-10;     %cm/ms
p.RTF = 25.690;
p.beta_pH = 0.04;
p.cap_0 = 1e-6;       %F/cm2
p.pH_C = 7.2;
p.psi_in = 0;          %mV
p.psi_out = -50;       %mV
p.q = 2.2;
p.r = 0.1;
p.tau_act = 1000;     %ms
p.tau_deact = 250;    %ms
p.Cl_C = 0.01;        %M
p.B = 0.0809542751;

p.init_V = 1.64636595e-16*N_Lys;    %L
p.S = 1.45267584e-8*N_Lys;          %cm2
p.N_CAX = 300*N_Lys;                   %number of CAX per lysosome
p.N_VATP = 300*N_Lys;                  %number of V-ATPase per lysosome

%Made up values. Can be modified
p.n_TPC = 1;           %Hill Coefficient
p.v_TPC = 1e-5 * 50 * N_Lys;      %M/ms
p.Ka_TPC = 0.03;      %uM, so about 30 nM

p.v_H = 6e-8*1e3;         %1/ms
p.v_Ca = p.v2; %1.49e-8*1e3;  %1/ms
p.psi_cax = -40;



%V-ATPase velocity
load('v_flux.mat');
p.v_flux = v_flux; %proton/s

save('parameters.mat', "p");
endfunction
