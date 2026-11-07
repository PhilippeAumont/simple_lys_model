function dS = Misc_channels(C_cav, V, Na_i, Ca_i, E_Na, f_cav_PLM_p, p)
  %Na-K pump and Phospholemman (PLM)
  % Parameters
  I_max_NaK = 4.0;           %pA/pF
  Km_Nai_np = 18800;        %uM
  Km_Nai_p = 13600;         %uM
  Km_Ko = 1500;             %uM
  k_PLM_PKA = 3.053e-6;    %1/uM ms
  K_PLM_PKA = 0.0011001;
  k_PLM_PP = 1.8491e-5;   %1/uM ms
  K_PLM_PP = 5.7392;

  PP_cav = 0.2;

  % Calculations
  df_cav_PLMp = k_PLM_PKA*C_cav*(1-f_cav_PLM_p)/(K_PLM_PKA+(1-f_cav_PLM_p))...
                - k_PLM_PP*PP_cav*f_cav_PLM_p/(K_PLM_PP+f_cav_PLM_p);

  sig = (e^(p.Na_o/67300) - 1)/7;
  f_NaK = inv(1 + 0.1245*e^(-0.1*V*p.F/(p.R*p.T)) + 0.0365*sig*e^(-V*p.F/(p.R*p.T)));
  Km_Nai = Km_Nai_np*(1-f_cav_PLM_p) + Km_Nai_p*f_cav_PLM_p;

  I_NaK = I_max_NaK*f_NaK*inv(1+(Km_Nai/Na_i)^3)*p.K_o/(p.K_o+Km_Ko);

  %Calcium pump
  I_pCa = p.I_max_pCa*Ca_i^2/(p.Km_pCa^2 + Ca_i^2);

  %Na/Ca Exchanger (NCX)
  I_NaCa = p.k_NaCa*inv(p.Km_Na^3 + p.Na_o^3)*inv(p.Km_Ca + p.Ca_o)...
          *inv(1 + p.k_sat*e^((p.n-1)*V*p.F/(p.R*p.T)))*((e^(p.n*V*p.F/(p.R*p.T))*Na_i^3*p.Ca_o)...
          - 2.0*e^((p.n-1)*V*p.F/(p.R*p.T))*p.Na_o^3*Ca_i);

  %Calcium Background
  E_CaN = (p.R*p.T/(2*p.F)) * log(p.Ca_o/Ca_i);
  I_Cab = p.G_Cab*(V - E_CaN);

  %Sodium Background
  I_Nab = p.G_Nab*(V - E_Na);

  %Ca-activated Cl current
  O_ClCa = 0.2/(1 + e^(-(V-46.7)/7.8));
  I_ClCa = p.G_ClCa*O_ClCa*(Ca_i/(Ca_i + p.Km_Cl))*(V-p.E_Cl);

  %Package output
  dS = [df_cav_PLMp; I_NaK; I_pCa; I_NaCa; I_Cab; I_Nab; I_ClCa];

endfunction
