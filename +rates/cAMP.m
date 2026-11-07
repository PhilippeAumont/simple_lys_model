function dcAMP = cAMP(S, cAMP, p);
  %Unpack state vectors
  dcAMP_cav_AC56 = S(1);
  dcAMP_ecav_AC47 = S(2);
  dcAMP_cyt_AC56 = S(3);
  dcAMP_cyt_AC47 = S(4);
  dcAMP_cav_PDE2 = S(5);
  dcAMP_cav_PDE3 = S(6);
  dcAMP_cav_PDE4 = S(7);
  dcAMP_ecav_PDE2 = S(8);
  dcAMP_ecav_PDE4 = S(9);
  dcAMP_cyt_PDE2 = S(10);
  dcAMP_cyt_PDE3 = S(11);
  dcAMP_cyt_PDE4 = S(12);
  dcAMP_cav_PKA = S(13);
  dcAMP_ecav_PKA = S(14);
  dcAMP_cyt_PKA = S(15);

  cAMP_cav = cAMP(1);
  cAMP_ecav = cAMP(2);
  cAMP_cyt = cAMP(3);

  %Parameters
  J_cav_ecav = 5.000e-12;  %uL/ms
  J_cav_cyt = 7.500e-11;   %uL/ms
  J_ecav_cyt = 9.000e-12;  %uL/ms

  %ODEs compilations
  dcAMP_cav = dcAMP_cav_PKA + dcAMP_cav_AC56 - dcAMP_cav_PDE2 - dcAMP_cav_PDE3...
             - dcAMP_cav_PDE4 - J_cav_ecav*(cAMP_cav-cAMP_ecav)/p.V_cav...
             - J_cav_cyt*(cAMP_cav-cAMP_cyt)/p.V_cav;
  dcAMP_ecav = dcAMP_ecav_PKA + dcAMP_ecav_AC47 - dcAMP_ecav_PDE2...
             - dcAMP_ecav_PDE4 - J_cav_ecav*(cAMP_ecav-cAMP_cav)/p.V_ecav...
             - J_ecav_cyt*(cAMP_ecav-cAMP_cyt)/p.V_ecav;
  dcAMP_cyt = dcAMP_cyt_PKA + dcAMP_cyt_AC56 + dcAMP_cyt_AC47- dcAMP_cyt_PDE2 - dcAMP_cyt_PDE3...
             - dcAMP_cyt_PDE4 - J_cav_cyt*(cAMP_cyt-cAMP_cav)/p.V_cyt...
             - J_ecav_cyt*(cAMP_cyt-cAMP_ecav)/p.V_cyt;


  dcAMP = [dcAMP_cav; dcAMP_ecav; dcAMP_cyt];
endfunction
