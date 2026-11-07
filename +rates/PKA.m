function dPKA = PKA(S, C, cAMP, p)
  %Unpack state vectors
  ARC_cav = S(1);
  A2RC_cav = S(2);
  A2R_cav = S(3);
  C_cav = S(4);
  PKIC_cav = S(5);
  ARC_ecav = S(6);
  A2RC_ecav = S(7);
  A2R_ecav = S(8);
  C_ecav = S(9);
  PKIC_ecav = S(10);
  ARC_cyt = S(11);
  A2RC_cyt = S(12);
  A2R_cyt = S(13);
  C_cyt = S(14);
  PKIC_cyt = S(15);

  C_cav = C(1);
  C_ecav = C(2);
  C_cyt = C(3);

  cAMP_cav = cAMP(1);
  cAMP_ecav = cAMP(2);
  cAMP_cyt = cAMP(3);

  %Parameters
  PKA_tot = 0.5176;                     %uM
  f_cav_PKA = 0.08;
  f_ecav_PKA = 0.20;
  f_cyt_PKA = 1-f_cav_PKA-f_ecav_PKA;
  PKI_tot = 0.4*PKA_tot;                 %uM
  f_cav_PKI = f_cav_PKA;
  f_ecav_PKI = f_ecav_PKA;
  f_cyt_PKI = f_cyt_PKA;
  k_PKAI_f1 = 5.6e-3;                  %1/uM ms
  K_PKAI_1 = 2.9;                       %uM
  k_PKAI_f2 = k_PKAI_f1;                 %1/uM ms
  K_PKAI_2 = 2.9;                       %uM
  k_PKAI_f3 = 2.6e-3;                  %1/ms
  K_PKAI_3 = 1.3;                       %uM
  k_PKI_f = 0.050;                     %1/uM ms
  K_PKI = 2.6e-4;                      %uM
  k_PKAII_f1 = k_PKAI_f1;                %1/uM ms
  K_PKAII_1 = 2.5;                      %uM
  k_PKAII_f2 = k_PKAI_f1;                %1/uM ms
  K_PKAII_2 = 2.5;                      %uM
  k_PKAII_f3 = k_PKAI_f3;                %1/ms
  K_PKAII_3 = K_PKAI_3;                  %uM

  %Concentrations and calculations
  %Caveolae (PKAI isoform)
  PKA_cav = f_cav_PKA*PKA_tot*p.V_cell/p.V_cav;
  RC_cav = 2*PKA_cav - ARC_cav - A2RC_cav - A2R_cav;
  PKI_cav = f_cav_PKI*PKI_tot*p.V_cell/p.V_cav - PKIC_cav;
  k_PKAII_b1 = k_PKAII_f1*K_PKAII_1;
  k_PKAII_b2 = k_PKAII_f2*K_PKAII_2;
  k_PKAII_b3 = k_PKAII_f3/K_PKAII_3;
  k_PKI_b = k_PKI_f*K_PKI;

  %Extracaveolae (PKAI isoform)
  PKA_ecav = f_ecav_PKA*PKA_tot*p.V_cell/p.V_ecav;
  RC_ecav = 2*PKA_ecav - ARC_ecav - A2RC_ecav - A2R_ecav;
  PKI_ecav = f_ecav_PKI*PKI_tot*p.V_cell/p.V_ecav - PKIC_ecav;

  %Cytosol  (PKAII isoform)
  PKA_cyt = f_cyt_PKA*PKA_tot*p.V_cell/p.V_cyt;
  RC_cyt = 2*PKA_cyt - ARC_cyt - A2RC_cyt - A2R_cyt;
  PKI_cyt = f_cyt_PKI*PKI_tot*p.V_cell/p.V_cyt - PKIC_cyt;
  k_PKAI_b1 = k_PKAI_f1*K_PKAI_1;
  k_PKAI_b2 = k_PKAI_f2*K_PKAI_2;
  k_PKAI_b3 = k_PKAI_f3/K_PKAI_3;

  %ODEs
  %Caveolae
  dcAMP_cav_PKA = -k_PKAII_f1*RC_cav*cAMP_cav + k_PKAII_b1*ARC_cav - k_PKAII_f2*ARC_cav*cAMP_cav...
                  + k_PKAII_b2*A2RC_cav;
  dARC_cav = k_PKAII_f1*RC_cav*cAMP_cav - k_PKAII_b1*ARC_cav - k_PKAII_f2*ARC_cav*cAMP_cav...
              + k_PKAII_b2*A2RC_cav;
  dA2RC_cav = k_PKAII_f2*ARC_cav*cAMP_cav - (k_PKAII_b2+k_PKAII_f3)*A2RC_cav...
              + k_PKAII_b3*A2R_cav*C_cav;
  dA2R_cav = k_PKAII_f3*A2RC_cav - k_PKAII_b3*A2R_cav*C_cav;
  dC_cav = k_PKAII_f3*A2RC_cav - k_PKAII_b3*A2R_cav*C_cav + k_PKI_b*PKIC_cav...
           - k_PKI_f*PKI_cav*C_cav;
  dPKIC_cav = -k_PKI_b*PKIC_cav + k_PKI_f*PKI_cav*C_cav;

  %Extracaveolae
  dcAMP_ecav_PKA = -k_PKAII_f1*RC_ecav*cAMP_ecav + k_PKAII_b1*ARC_ecav - k_PKAII_f2*ARC_ecav*cAMP_ecav...
                  + k_PKAII_b2*A2RC_ecav;
  dARC_ecav = k_PKAII_f1*RC_ecav*cAMP_ecav - k_PKAII_b1*ARC_ecav - k_PKAII_f2*ARC_ecav*cAMP_ecav...
              + k_PKAII_b2*A2RC_ecav;
  dA2RC_ecav = k_PKAII_f2*ARC_ecav*cAMP_ecav - (k_PKAII_b2+k_PKAII_f3)*A2RC_ecav...
              + k_PKAII_b3*A2R_ecav*C_ecav;
  dA2R_ecav = k_PKAII_f3*A2RC_ecav - k_PKAII_b3*A2R_ecav*C_ecav;
  dC_ecav = k_PKAII_f3*A2RC_ecav - k_PKAII_b3*A2R_ecav*C_ecav + k_PKI_b*PKIC_ecav...
           - k_PKI_f*PKI_ecav*C_ecav;
  dPKIC_ecav = -k_PKI_b*PKIC_ecav + k_PKI_f*PKI_ecav*C_ecav;

  %Cytosol
  dcAMP_cyt_PKA = -k_PKAI_f1*RC_cyt*cAMP_cyt + k_PKAI_b1*ARC_cyt - k_PKAI_f2*ARC_cyt*cAMP_cyt...
                  + k_PKAI_b2*A2RC_cyt;
  dARC_cyt = k_PKAI_f1*RC_cyt*cAMP_cyt - k_PKAI_b1*ARC_cyt - k_PKAI_f2*ARC_cyt*cAMP_cyt...
              + k_PKAI_b2*A2RC_cyt;
  dA2RC_cyt = k_PKAI_f2*ARC_cyt*cAMP_cyt - (k_PKAI_b2+k_PKAI_f3)*A2RC_cyt...
              + k_PKAI_b3*A2R_cyt*C_cyt;
  dA2R_cyt = k_PKAI_f3*A2RC_cyt - k_PKAI_b3*A2R_cyt*C_cyt;
  dC_cyt = k_PKAI_f3*A2RC_cyt - k_PKAI_b3*A2R_cyt*C_cyt + k_PKI_b*PKIC_cyt...
           - k_PKI_f*PKI_cyt*C_cyt;
  dPKIC_cyt = -k_PKI_b*PKIC_cyt + k_PKI_f*PKI_cyt*C_cyt;

  %Package output
  dPKA = [dcAMP_cav_PKA;dARC_cav;dA2RC_cav;dA2R_cav;dC_cav;dPKIC_cav;dcAMP_ecav_PKA;
          dARC_ecav;dA2RC_ecav;dA2R_ecav;dC_ecav;dPKIC_ecav;dcAMP_cyt_PKA;dARC_cyt;
          dA2RC_cyt;dA2R_cyt;dC_cyt;dPKIC_cyt];

endfunction
