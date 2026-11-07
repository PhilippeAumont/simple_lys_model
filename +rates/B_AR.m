function dB_AR = B_AR(S, C, p)

  %Unpack state vector
  R_cav_PKA = S(1);
  R_cav_GRK2 = S(2);
  Gs_cav_aGTP = S(3);
  Gs_cav_By = S(4);
  Gs_cav_aGDP = S(5);
  R_ecav_PKA = S(6);
  R_ecav_GRK2 = S(7);
  Gs_ecav_aGTP = S(8);
  Gs_ecav_By = S(9);
  Gs_ecav_aGDP = S(10);
  R_cyt_PKA = S(11);
  R_cyt_GRK2 = S(12);
  Gs_cyt_aGTP = S(13);
  Gs_cyt_By = S(14);
  Gs_cyt_aGDP = S(15);

  C_cav = C(1);
  C_ecav = C(2);
  C_cyt = C(3);

  %Parameters
  R_tot = 0.0103;                    %uM
  f_cav_B1 = 0.01;
  f_ecav_B1 = 0.5;
  f_cyt_B1 = 1-f_cav_B1-f_ecav_B1;
  Gs_tot = 2.054;                    %uM
  f_cav_Gs = 0.4;
  f_ecav_Gs = 0.4;
  f_cyt_Gs = 1-f_cav_Gs-f_ecav_Gs;
  K_B1_L = 0.567;                    %uM
  K_B1_H = 0.0617;                   %uM
  K_B1_C = 2.86;                     %uM
  k_PKA_on = 8.1e-07;               %1/uM ms
  k_PKA_off = 2.025e-07;             %1/ms
  k_GRK2_on = 2.43e-07;             %1/ms
  k_GRK2_off = k_PKA_off;               %1/ms
  k_act1_Gs = 4.9e-3;                %1/ms
  k_act2_Gs = 2.6e-4;               %1/ms
  k_hyd_Gs = 8.0e-04;               %1/ms
  k_reas_Gs = 1.2;                   %1/uM ms

  %Concentrations
  %Caveolae (cav)
  R_cav_tot = f_cav_B1*R_tot*p.V_cell/p.V_cav;
  Gs_cav_aBy = f_cav_Gs*Gs_tot*p.V_cell/p.V_cav - Gs_cav_aGTP - Gs_cav_aGDP;
  R_cav_np = R_cav_tot - R_cav_PKA - R_cav_GRK2;
  a_cav_B1 = inv(K_B1_L)*(K_B1_L+p.L)*(K_B1_H+p.L);
  b_cav_B1 = Gs_cav_aBy*(K_B1_H+p.L)-(R_cav_np*(K_B1_H+p.L))+(K_B1_C*K_B1_H*(1+p.L/K_B1_L));
  c_cav_B1 = -R_cav_np*K_B1_C*K_B1_H;
  R_cav_np_f = (-b_cav_B1+sqrt(b_cav_B1^2 - 4*a_cav_B1*c_cav_B1))/(2*a_cav_B1);
  Gs_cav_f = Gs_cav_aBy*inv(1+R_cav_np_f*(inv(K_B1_C) + p.L/(K_B1_C*K_B1_H)));
  LR_cav_np = p.L*R_cav_np_f/K_B1_L;
  R_Gs_cav_np = R_cav_np_f*Gs_cav_f/K_B1_C;
  LR_Gs_cav_np = p.L*R_cav_np_f*Gs_cav_f/(K_B1_C*K_B1_H);

  %Extracaveolae (ecav)
  R_ecav_tot = f_ecav_B1*R_tot*p.V_cell/p.V_ecav;
  Gs_ecav_aBy = f_ecav_Gs*Gs_tot*p.V_cell/p.V_ecav - Gs_ecav_aGTP - Gs_ecav_aGDP;
  R_ecav_np = R_ecav_tot - R_ecav_PKA - R_ecav_GRK2;
  a_ecav_B1 = a_cav_B1;
  b_ecav_B1 = Gs_ecav_aBy*(K_B1_H+p.L)-(R_ecav_np*(K_B1_H+p.L))+(K_B1_C*K_B1_H*(1+p.L/K_B1_L));
  c_ecav_B1 = -R_ecav_np*K_B1_C*K_B1_H;
  R_ecav_np_f = (-b_ecav_B1+sqrt(b_ecav_B1^2 - 4*a_ecav_B1*c_ecav_B1))/(2*a_ecav_B1);
  Gs_ecav_f = Gs_ecav_aBy*inv(1+R_ecav_np_f*(inv(K_B1_C) + p.L/(K_B1_C*K_B1_H)));
  LR_ecav_np = p.L*R_ecav_np_f/K_B1_L;
  R_Gs_ecav_np = R_ecav_np_f*Gs_ecav_f/K_B1_C;
  LR_Gs_ecav_np = p.L*R_ecav_np_f*Gs_ecav_f/(K_B1_C*K_B1_H);

  %Cytosol (cyt)
  R_cyt_tot = f_cyt_B1*R_tot*p.V_cell/p.V_cyt;
  Gs_cyt_aBy = f_cyt_Gs*Gs_tot*p.V_cell/p.V_cyt - Gs_cyt_aGTP - Gs_cyt_aGDP;
  R_cyt_np = R_cyt_tot - R_cyt_PKA - R_cyt_GRK2;
  a_cyt_B1 = a_cav_B1;
  b_cyt_B1 = Gs_cyt_aBy*(K_B1_H+p.L)-(R_cyt_np*(K_B1_H+p.L))+(K_B1_C*K_B1_H*(1+p.L/K_B1_L));
  c_cyt_B1 = -R_cyt_np*K_B1_C*K_B1_H;
  R_cyt_np_f = (-b_cyt_B1+sqrt(b_cyt_B1^2 - 4*a_cyt_B1*c_cyt_B1))/(2*a_cyt_B1);
  Gs_cyt_f = Gs_cyt_aBy*inv(1+R_cyt_np_f*(inv(K_B1_C) + p.L/(K_B1_C*K_B1_H)));
  LR_cyt_np = p.L*R_cyt_np_f/K_B1_L;
  R_Gs_cyt_np = R_cyt_np_f*Gs_cyt_f/K_B1_C;
  LR_Gs_cyt_np = p.L*R_cyt_np_f*Gs_cyt_f/(K_B1_C*K_B1_H);


  %dXdt
  dR_cav_PKA = k_PKA_on*C_cav*R_cav_np - k_PKA_off*R_cav_PKA;
  dR_cav_GRK2 = k_GRK2_on*(LR_cav_np + LR_Gs_cav_np) - k_GRK2_off*R_cav_GRK2;
  dGs_cav_aGTP = k_act2_Gs*R_Gs_cav_np + k_act1_Gs*LR_Gs_cav_np - k_hyd_Gs*Gs_cav_aGTP;
  dGs_cav_By = k_act2_Gs*R_Gs_cav_np + k_act1_Gs*LR_Gs_cav_np - k_reas_Gs*Gs_cav_By*Gs_cav_aGDP;
  dGs_cav_aGDP = k_hyd_Gs*Gs_cav_aGTP - k_reas_Gs*Gs_cav_By*Gs_cav_aGDP;

  dR_ecav_PKA = k_PKA_on*C_ecav*R_ecav_np - k_PKA_off*R_ecav_PKA;
  dR_ecav_GRK2 = k_GRK2_on*(LR_ecav_np + LR_Gs_ecav_np) - k_GRK2_off*R_ecav_GRK2;
  dGs_ecav_aGTP = k_act2_Gs*R_Gs_ecav_np + k_act1_Gs*LR_Gs_ecav_np - k_hyd_Gs*Gs_ecav_aGTP;
  dGs_ecav_By = k_act2_Gs*R_Gs_ecav_np + k_act1_Gs*LR_Gs_ecav_np - k_reas_Gs*Gs_ecav_By*Gs_ecav_aGDP;
  dGs_ecav_aGDP = k_hyd_Gs*Gs_ecav_aGTP - k_reas_Gs*Gs_ecav_By*Gs_ecav_aGDP;

  dR_cyt_PKA = k_PKA_on*C_cyt*R_cyt_np - k_PKA_off*R_cyt_PKA;
  dR_cyt_GRK2 = k_GRK2_on*(LR_cyt_np + LR_Gs_cyt_np) - k_GRK2_off*R_cyt_GRK2;
  dGs_cyt_aGTP = k_act2_Gs*R_Gs_cyt_np + k_act1_Gs*LR_Gs_cyt_np - k_hyd_Gs*Gs_cyt_aGTP;
  dGs_cyt_By = k_act2_Gs*R_Gs_cyt_np + k_act1_Gs*LR_Gs_cyt_np - k_reas_Gs*Gs_cyt_By*Gs_cyt_aGDP;
  dGs_cyt_aGDP = k_hyd_Gs*Gs_cyt_aGTP - k_reas_Gs*Gs_cyt_By*Gs_cyt_aGDP;

  dB_AR = [dR_cav_PKA; dR_cav_GRK2; dGs_cav_aGTP; dGs_cav_By; dGs_cav_aGDP;
             dR_ecav_PKA;dR_ecav_GRK2;dGs_ecav_aGTP;dGs_ecav_By;dGs_ecav_aGDP;
             dR_cyt_PKA; dR_cyt_GRK2; dGs_cyt_aGTP; dGs_cyt_By; dGs_cyt_aGDP];

 endfunction
