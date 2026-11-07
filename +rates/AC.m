function dAC = AC(S, p)
  %Unpack state vectorize
  Gs_cav_aGTP = S(1);
  Gs_cav_By = S(2);
  Gs_ecav_aGTP = S(3);
  Gs_ecav_By = S(4);
  Gs_cyt_aGTP = S(5);
  Gs_cyt_By = S(6);

  %Parameters
  K_m_ATP = 340;            %uM
  ATP = 5000;               %uM
  AC_tot = 0.02622;        %uM
  f_56 = 0.74;
  f_cav_56 = 0.0875;
  f_ecav_47 = 0.1648;
  K_m_56_Gsa = 0.0852;      %uM
  h_56_Gsa = 1.357;
  V_56_GBy = 1.430;
  K_m_56_GsBy = 0.003793;  %uM
  h_56_GsBy = 1.0842;
  AC56_basal = 0.0377;
  AF_56 = 0.051133;        %1/ms
  K_m_47_Gsa = 0.05008;     %uM
  h_47_Gsa = 1.1657;
  V_47_GBy = 1.3500;
  K_m_47_GsBy = 0.004466;  %uM
  h_47_GsBy = 0.8700;
  AC47_basal = 0.04725;
  AF_47 = 9.2830e-03;     %1/s

  %Concentrations and calculations
  AC56_cav = f_cav_56*f_56*AC_tot*p.V_cell/p.V_cav;
  k_cav_56 = AF_56*(AC56_basal + Gs_cav_aGTP^h_56_Gsa/(K_m_56_Gsa+Gs_cav_aGTP^h_56_Gsa))...
            *(1 + V_56_GBy*Gs_cav_By^h_56_GsBy/(K_m_56_GsBy + Gs_cav_By^h_56_GsBy));

  AC47_ecav = f_ecav_47*(1-f_56)*AC_tot*p.V_cell/p.V_ecav;
  k_ecav_47 = AF_47*(AC47_basal + Gs_ecav_aGTP^h_47_Gsa/(K_m_47_Gsa+Gs_ecav_aGTP^h_47_Gsa))...
            *(1 + V_47_GBy*Gs_ecav_By^h_47_GsBy/(K_m_47_GsBy + Gs_ecav_By^h_47_GsBy));

  AC56_cyt = (1-f_cav_56)*f_56*AC_tot*p.V_cell/p.V_cyt;
  AC47_cyt = (1-f_ecav_47)*(1-f_56)*AC_tot*p.V_cell/p.V_cyt;
  k_cyt_56 = AF_56*(AC56_basal + Gs_cyt_aGTP^h_56_Gsa/(K_m_56_Gsa+Gs_cyt_aGTP^h_56_Gsa))...
            *(1 + V_56_GBy*Gs_cyt_By^h_56_GsBy/(K_m_56_GsBy + Gs_cyt_By^h_56_GsBy));
  k_cyt_47 = AF_47*(AC47_basal + Gs_cyt_aGTP^h_47_Gsa/(K_m_47_Gsa+Gs_cyt_aGTP^h_47_Gsa))...
            *(1 + V_47_GBy*Gs_cyt_By^h_47_GsBy/(K_m_47_GsBy + Gs_cyt_By^h_47_GsBy));

  %dXdt
  dcAMP_cav_56 = k_cav_56*AC56_cav*ATP/(K_m_ATP + ATP);
  dcAMP_ecav_47 = k_ecav_47*AC47_ecav*ATP/(K_m_ATP + ATP);
  dcAMP_cyt_56 = k_cyt_56*AC56_cyt*ATP/(K_m_ATP + ATP);
  dcAMP_cyt_47 = k_cyt_47*AC47_cyt*ATP/(K_m_ATP + ATP);

  dAC = [dcAMP_cav_56;dcAMP_ecav_47;dcAMP_cyt_56;dcAMP_cyt_47];
endfunction
