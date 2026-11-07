function dC = Concentrations(Ca, J, I, p)
  %Unpack state vector
  Ca_i = Ca(1);
  Ca_ss = Ca(2);
  Ca_md = Ca(3);
  Ca_JSR = Ca(4);

  J_rel = J(1);
  J_rel_md = J(2);
  J_tr = J(3);
  J_xfer = J(4);
  J_xfer_md = J(5);
  J_leak = J(6);
  J_up = J(7);
  J_trpn = J(8);
  J_Ca = J(9)*(p.V_lys/p.V_cyt);
  J_CAX = J(10)*(p.V_lys/p.V_cyt);
  J_TPC = J(11)*(p.V_lys/p.V_md);

  I_cav_CaL = I(1);
  I_ecav_CaL = I(2);
  I_Na = I(3);
  I_Kr = I(4);
  I_Kur = I(5);
  I_Kto_f = I(6);
  I_K1 = I(7);
  I_Kss = I(8);
  I_NaK = I(9);
  I_pCa = I(10);
  I_NaCa = I(11);
  I_Cab = I(12);
  I_Nab = I(13);
  I_ClCa = I(14);

  I_Kto_s = 0;  %These are not in the model. Likely they are assumed to be 0.
  I_Ks = 0;

  %Calcium
  B_i = 1/(1 + p.CMDN_tot*p.Km_CMDN/(p.Km_CMDN+Ca_i)^2);
  B_ss = 1/(1 + p.CMDN_tot*p.Km_CMDN/(p.Km_CMDN+Ca_ss)^2);
  B_JSR = 1/(1 + p.CSQN_tot*p.Km_CSQN/(p.Km_CSQN+Ca_JSR)^2);
  B_md = 1/(1 + p.CMDN_tot*p.Km_CMDN/(p.Km_CMDN+Ca_md)^2);

  dCa_i = B_i*(J_leak + J_xfer + J_xfer_md - J_up - J_trpn - J_Ca - J_CAX - (I_Cab - 2*I_NaCa + I_pCa...
                + I_cav_CaL)*(p.A_cap*p.C_m/(2*p.V_cyt*p.F)));
  dCa_ss = B_ss*(J_rel*p.V_JSR/p.V_ss - J_xfer*p.V_cyt/p.V_ss...
          - I_ecav_CaL*p.A_cap*p.C_m/(2*p.V_ss*p.F));
  dCa_JSR = B_JSR*(J_tr - J_rel);
  dCa_NSR = (J_up - J_leak)*p.V_cyt/p.V_NSR - J_tr*p.V_JSR/p.V_NSR;
  dCa_md = B_md*(J_rel_md*p.V_NSR/p.V_md - J_TPC - J_xfer_md*p.V_cyt/p.V_md);

  %Sodium
  dNa_i = -(I_Na + I_Nab + 3*I_NaCa + 3*I_NaK)*p.A_cap*p.C_m/(p.V_cyt*p.F);

  %Potassium
  dK_i = -(I_Kto_f + I_Kto_s + I_Kur + I_Kss + I_K1 + I_Kr + I_Ks - 2*I_NaK)...
            *p.A_cap*p.C_m/(p.V_cyt*p.F);

  %Pack outputs
  dC = [dCa_i; dCa_ss; dCa_md; dCa_JSR; dCa_NSR; dNa_i; dK_i];

endfunction
