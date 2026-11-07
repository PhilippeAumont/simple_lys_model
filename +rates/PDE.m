function dPDE = PDE(S, C, cAMP, p)
  %Unpack state vector
  PDE3_cav_p = S(1);
  PDE4_cav_p = S(2);
  PDE4_ecav_p = S(3);
  PDE3_cyt_p = S(4);
  PDE4_cyt_p = S(5);

  C_cav = C(1);
  C_ecav = C(2);
  C_cyt = C(3);

  cAMP_cav = cAMP(1);
  cAMP_ecav = cAMP(2);
  cAMP_cyt = cAMP(3);

  %Parameters
  h_IBMX_PDE2 = 1.000;
  K_IBMX_PDE2 = 29.50;    %uM
  h_IBMX_PDE3 = 1.000;
  K_IBMX_PDE3 = 5.100;    %uM
  h_IBMX_PDE4 = 1.000;
  K_IBMX_PDE4 = 16.200;   %uM
  kf_PDEp = 1.96e-05;     %1/uM ms
  kb_PDEp = 1.02e-05;     %1/ms
  d_k_PDE34 = 3.0;
  k_PDE2 = 0.020;         %1/ms
  K_m_PDE2 = 33;           %uM
  k_PDE3 = 2.5e-3;       %1/ms
  K_m_PDE3 = 0.44;        %uM
  k_PDE4 = 3.5e-3;       %1/ms
  K_m_PDE4 = 1.4;         %uM
  f_PDE_part = 0.2;
  r_part_PDE23 = 0.570;
  r_part_PDE34 = 0.748;
  PDE2_tot = 0.034610;    %uM
  PDE3_tot = 0.010346;    %uM
  PDE4_tot = 0.026687;    %uM
  f_cav_PDE2 = 0.06608;
  f_ecav_PDE2 = 2*f_cav_PDE2;
  f_cyt_PDE2 = 1-f_cav_PDE2-f_ecav_PDE2;
  f_cav_PDE3 = 0.29814;
  f_ecav_PDE3 = 0.0;
  f_cyt_PDE3 = 1-f_cav_PDE3-f_ecav_PDE3;
  f_cav_PDE4 = 0.05366;
  f_ecav_PDE4 = 2*f_cav_PDE4;
  f_cyt_PDE4 = 1-f_cav_PDE4-f_ecav_PDE4;

  %Concentrations and calculations
  %Caveolae (cav)
  PDE2_cav = (1 - p.IBMX^h_IBMX_PDE2/(K_IBMX_PDE2+p.IBMX^h_IBMX_PDE2))*f_cav_PDE2*PDE2_tot*...
              p.V_cell/p.V_cav;
  PDE3_cav = (1 - p.IBMX^h_IBMX_PDE3/(K_IBMX_PDE3+p.IBMX^h_IBMX_PDE3))*f_cav_PDE3*PDE3_tot*...
              p.V_cell/p.V_cav;
  PDE4_cav = (1 - p.IBMX^h_IBMX_PDE4/(K_IBMX_PDE4+p.IBMX^h_IBMX_PDE4))*f_cav_PDE4*PDE4_tot*...
              p.V_cell/p.V_cav;

  %Extracaveolae (ecav)
  PDE2_ecav = (1 - p.IBMX^h_IBMX_PDE2/(K_IBMX_PDE2+p.IBMX^h_IBMX_PDE2))*f_ecav_PDE2*PDE2_tot*...
              p.V_cell/p.V_ecav;
  PDE4_ecav = (1 - p.IBMX^h_IBMX_PDE4/(K_IBMX_PDE4+p.IBMX^h_IBMX_PDE4))*f_ecav_PDE4*PDE4_tot*...
              p.V_cell/p.V_ecav;

  %Cytosol (cyt)
  PDE2_cyt = (1 - p.IBMX^h_IBMX_PDE2/(K_IBMX_PDE2+p.IBMX^h_IBMX_PDE2))*f_cyt_PDE2*PDE2_tot*...
              p.V_cell/p.V_cyt;
  PDE3_cyt = (1 - p.IBMX^h_IBMX_PDE3/(K_IBMX_PDE3+p.IBMX^h_IBMX_PDE3))*f_cyt_PDE3*PDE3_tot*...
              p.V_cell/p.V_cyt;
  PDE4_cyt = (1 - p.IBMX^h_IBMX_PDE4/(K_IBMX_PDE4+p.IBMX^h_IBMX_PDE4))*f_cyt_PDE4*PDE4_tot*...
              p.V_cell/p.V_cyt;

  %ODE
  %Caveolae (cav)
  dPDE3_cav_p = kf_PDEp*C_cav*(PDE3_cav-PDE3_cav_p) - kb_PDEp*PDE3_cav_p;
  dPDE4_cav_p = kf_PDEp*C_cav*(PDE4_cav-PDE4_cav_p) - kb_PDEp*PDE4_cav_p;
  dcAMP_cav_PDE2 = k_PDE2*PDE2_cav*cAMP_cav/(K_m_PDE2+cAMP_cav);
  dcAMP_cav_PDE3 = (k_PDE3*(PDE3_cav-PDE3_cav_p)*cAMP_cav + d_k_PDE34*k_PDE3*PDE3_cav_p*cAMP_cav)...
                    /(K_m_PDE3+cAMP_cav);
  dcAMP_cav_PDE4 = (k_PDE4*(PDE4_cav-PDE4_cav_p)*cAMP_cav + d_k_PDE34*k_PDE4*PDE4_cav_p*cAMP_cav)...
                    /(K_m_PDE4+cAMP_cav);

  %Extracaveolae (ecav)
  dPDE4_ecav_p = kf_PDEp*C_ecav*(PDE4_ecav-PDE4_ecav_p) - kb_PDEp*PDE4_ecav_p;
  dcAMP_ecav_PDE2 = k_PDE2*PDE2_ecav*cAMP_ecav/(K_m_PDE2+cAMP_ecav);
  dcAMP_ecav_PDE4 = (k_PDE4*(PDE4_ecav-PDE4_ecav_p)*cAMP_ecav + d_k_PDE34*k_PDE4*PDE4_ecav_p*cAMP_ecav)...
                    /(K_m_PDE4+cAMP_ecav);

  %Cytosol (cyt)
  dPDE3_cyt_p = kf_PDEp*C_cyt*(PDE3_cyt-PDE3_cyt_p) - kb_PDEp*PDE3_cyt_p;
  dPDE4_cyt_p = kf_PDEp*C_cyt*(PDE4_cyt-PDE4_cyt_p) - kb_PDEp*PDE4_cyt_p;
  dcAMP_cyt_PDE2 = k_PDE2*PDE2_cyt*cAMP_cyt/(K_m_PDE2+cAMP_cyt);
  dcAMP_cyt_PDE3 = (k_PDE3*(PDE3_cyt-PDE3_cyt_p)*cAMP_cyt + d_k_PDE34*k_PDE3*PDE3_cyt_p*cAMP_cyt)...
                    /(K_m_PDE3+cAMP_cyt);
  dcAMP_cyt_PDE4 = (k_PDE4*(PDE4_cyt-PDE4_cyt_p)*cAMP_cyt + d_k_PDE34*k_PDE4*PDE4_cyt_p*cAMP_cyt)...
                    /(K_m_PDE4+cAMP_cyt);

  dPDE = [dPDE3_cav_p;dPDE4_cav_p;dcAMP_cav_PDE2;dcAMP_cav_PDE3;dcAMP_cav_PDE4;
          dPDE4_ecav_p;dcAMP_ecav_PDE2;dcAMP_ecav_PDE4;dPDE3_cyt_p;dPDE4_cyt_p;
          dcAMP_cyt_PDE2;dcAMP_cyt_PDE3;dcAMP_cyt_PDE4];
endfunction
