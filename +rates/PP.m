function dPP = PP(S, C_cyt)
  %Unpack state vectorize
  Inhib1_cyt_p_tot = S(1);

  %Parameters
  PP1_cyt = 0.2;           %uM
  PP2A_cyt = 0.0607843;  %uM
  PP1_cav = 0.1;           %uM
  PP2A_cav = 0.1;          %uM
  PP_cav = 0.2;            %uM
  PP1_ecav = 0.1;          %uM
  Inhib1_cyt = 0.08543;   %uM
  K_inh1 = 1.0e-3;        %uM
  k_PKA_inh1 = 1.0800;    %1/uM ms
  Km_PKA_inh1 = 1.5;       %uM
  k_PP2A_inh1 = 0.050670; %1/uM ms
  Km_PP2A_inh1 = 1.0e-3;  %uM

  %Concentrations and calculations
  Inhib1_cyt_f = Inhib1_cyt - Inhib1_cyt_p_tot;
  a_inh1 = 1.0;
  b_inh1 = K_inh1 + PP1_cyt - Inhib1_cyt_p_tot;
  c_inh1 = -Inhib1_cyt_p_tot*K_inh1;

  Inhib1_cyt_p = (-b_inh1 + sqrt(b_inh1^2 - 4*a_inh1*c_inh1))/(2*a_inh1);
  PP1_cyt_f = PP1_cyt*K_inh1/(K_inh1 + Inhib1_cyt_p); %Used in PLB Module

  %ODEs
  dInhib1_cyt_p_tot = (k_PKA_inh1*C_cyt*Inhib1_cyt_f)/(Km_PKA_inh1+Inhib1_cyt_f)...
                     - (k_PP2A_inh1*PP2A_cyt*Inhib1_cyt_p_tot)/(Km_PP2A_inh1 + Inhib1_cyt_p_tot);

  %Package output
  dPP = [dInhib1_cyt_p_tot; PP1_cyt_f];
endfunction
