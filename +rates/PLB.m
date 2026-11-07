function dS = PLB(S, C_cyt)
  %Unpack state vector
  f_cyt_PLB_p = S(1);
  PP1_cyt_f = S(2);

  %Parameters
  Km_up_np = 0.41;             %uM
  Km_up_p = 0.31;              %uM
  k_PLB_PKA = 0.108917e-4;   %1/uM ms
  K_PLB_PKA = 4.90970e-4;
  k_PLB_PP1 = 4.41956e-5;    %1/uM ms
  K_PLB_PP1 = 1.69376e-2;

  %Calculations
  Km_up = Km_up_np*(1-f_cyt_PLB_p) + Km_up_p*f_cyt_PLB_p;
  df_cyt_PLB_p = k_PLB_PKA*C_cyt*(1-f_cyt_PLB_p)/(K_PLB_PKA+(1-f_cyt_PLB_p))...
                - k_PLB_PP1*PP1_cyt_f*f_cyt_PLB_p/(K_PLB_PP1+f_cyt_PLB_p);

  %Pack outputs
  dS = [Km_up; df_cyt_PLB_p];
endfunction
