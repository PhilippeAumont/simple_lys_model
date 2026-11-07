function dS = K_channels(S, C_ecav, V, K_i, p)
  %Unpack state vector
  f_ecav_IKur = S(1);
  a_ur = S(2);
  i_ur = S(3);
  a_urp = S(4);
  i_urp = S(5);

  f_ecav_IKto = S(6);
  a_to_f = S(7);
  i_to_f = S(8);
  a_to_fp = S(9);
  i_to_fp = S(10);

  a_Kss = S(11);

  %Ultra-rapidly activating delayed rectifier K channel - Kur
  % Parameters
  G_Kur = 0.3424;                %pA/pF
  G_Kurp = 0.53307;              %pA/pF
  k_IKur_PKA = 6.9537e-6;       %1/uM ms
  K_IKur_PKA = 0.138115;
  k_IKur_PP = 3.170e-5;         %1/uM ms
  K_IKur_PP = 0.23310;

  PP1_ecav = 0.1;

  % Calculations
  E_K = (p.R*p.T/p.F)*log(p.K_o/K_i);
  I_Kur = (G_Kur*a_ur*i_ur*f_ecav_IKur + G_Kurp*a_urp*i_urp*(1-f_ecav_IKur))*(V-E_K);

  df_ecav_IKur = k_IKur_PP*PP1_ecav*(1-f_ecav_IKur)/(K_IKur_PP+(1-f_ecav_IKur))...
                - k_IKur_PKA*C_ecav*f_ecav_IKur/(K_IKur_PKA+f_ecav_IKur);

  a_ss = inv(1 + e^(-(V+22.5)/7.7));
  i_ss = inv(1 + e^(-(V+45.2)/5.7));
  t_aur = 6.1/(e^(0.0629*(V+40.0)) + e^(-0.0629*(V+40.0))) + 2.058;
  t_iur = 1200.0 - 170.0/(1 + e^((V+45.2)/5.7));

  da_ur = (a_ss-a_ur)/t_aur;
  di_ur = (i_ss - i_ur)/t_iur;
  da_urp = (a_ss - a_urp)/t_aur;
  di_urp = (i_ss - i_urp)/t_iur;

  %Rapidly inactivating transient outward K channel - Kto_f
  % Parameters
  G_Kto_f = 0.3846;          %pA/pF
  G_Kto_fp = G_Kto_f;          %pA/pF
  k_IKto_PKA = 4.38983e-5;  %1/uM ms
  K_IKto_PKA = 0.27623;
  k_IKto_PP = 9.09678e-5;   %1/uM ms
  K_IKto_PP = 0.23310;

  %Calculations
  I_Kto_f = (G_Kto_f*a_to_f^3*i_to_f*(1-f_ecav_IKto) + G_Kto_fp*a_to_fp^3*i_to_fp*f_ecav_IKto)...
              * (V-E_K);

  df_ecav_IKto = k_IKto_PKA*C_ecav*(1-f_ecav_IKto)/(K_IKto_PKA+(1-f_ecav_IKto))...
                - k_IKto_PP*PP1_ecav*f_ecav_IKto/(K_IKto_PP+f_ecav_IKto);

  a_a = 0.18064*e^(0.03577*(V+33.0));
  B_a = 0.3956*e^(-0.06237*(V+33.0));
  a_i = 0.000152*e^(-(V+15.5)/7.0)/(0.067083*e^(-(V+35.5)/7.0) + 1);
  B_i = 0.00095*e^((V+35.5)/7.0)/(0.051335*e^((V+35.5)/7.0) + 1);
  a_ap = 0.18064*e^(0.03577*(V+17.0));
  B_ap = 0.3956*e^(-0.06237*(V+17.0));
  a_ip = 0.000152*e^(-(V+7.5)/7.0)/(0.067083*e^(-(V+27.5)/7.0) + 1);
  B_ip = 0.00095*e^((V+27.5)/7.0)/(0.051335*e^((V+27.5)/7.0) + 1);

  da_to_f = a_a*(1-a_to_f) - B_a*a_to_f;
  di_to_f = a_i*(1-i_to_f) - B_i*i_to_f;
  da_to_fp = a_ap*(1-a_to_fp) - B_ap*a_to_fp;
  di_to_fp = a_ip*(1-i_to_fp) - B_ip*i_to_fp;

  %Other Channels
  %Time-independent K channel - K1
  a_K1 = 1.02/(1 + exp(0.2385*(V - E_K - 59.215)));
  B_K1 = (0.8*exp(0.08032*(V-E_K+5.476)) + exp(0.06175*(V-E_K-594.31)))...
          / (1 + exp(-0.5143*(V-E_K+4.753)));

  I_K1 = 0.27*sqrt(p.K_o/5400)*(a_K1/(a_K1 + B_K1))*(V-E_K);

  %Noninactivating Steady-State Voltage activated K current - Kss
  I_Kss = p.G_Kss*a_Kss*(V-E_K);
  t_Kss = 1235.5/(e^(0.0862*(V+40.0)) + e^(-0.0862*(V+40.0))) + 13.17;
  da_Kss = (a_ss-a_Kss)/t_Kss;

  %Pack output
  dS = [
  I_Kur; I_Kto_f; I_K1; I_Kss;
  da_to_f; di_to_f; da_ur; di_ur; da_Kss; df_ecav_IKur;
  da_urp; di_urp; df_ecav_IKto; da_to_fp; di_to_fp
  ];

endfunction
