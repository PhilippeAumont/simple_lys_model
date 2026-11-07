function dS_RyR = RyR_md(S, Ca_md, C_cyt, PP1_cyt_f, p)
  %Unpack state vector
  o1 = S(1);
  o2 = S(2);
  c2 = S(3);
  o1p = S(4);
  o2p = S(5);
  c1p = S(6);
  c2p = S(7);

  c1 = 1-(o1+o2+c2+o1p+o2p+c1p+c2p);

  Y = [o1;o2;c1;c2;o1p;o2p;c1p;c2p];

  %Parameters
  RyR_tot = 0.1993;        %uM
  n = 4;
  m = 3;
  k_on_a = 6.075e-3;       %1/uM^4 ms
  k_off_a = 0.07125*5;        %1/s
  k_on_b = 4.05e-3;        %1/uM^3 ms
  k_off_b = 0.9650;         %1/ms
  k_on_c = 9.0e-3;         %1/ms
  k_off_c = 8e-4*5;            %1/ms
  k_on_ap = 5*k_on_a;        %1/uM^4 ms
  k_off_ap = 3*k_off_a;        %1/ms
  k_on_bp = 5*k_on_b;        %1/uM^3 ms
  k_off_bp = 3*k_off_b;        %1/ms
  k_on_cp = 50*k_on_c;       %1/ms
  k_off_cp = 30*k_off_c;       %1/ms
  f_RyR = 0.001;
  k_RyR_PKA = 5.775e-5;    %1/uM ms
  k_RyR_PP = 2.8875e-04;  %1/uM ms
  K_RyR_PKA = 0.5;           %uM
  K_RyR_PP = 0.05;          %uM

  PP1_cyt = 0.2;           %uM
  PP2_cyt = 0.060784;    %uM

  PP_cyt = PP1_cyt_f + PP2_cyt;

  %Calculations
  RyR_md = 0.01*RyR_tot*p.V_cell/p.V_md;

  %Functions for phosphorylation and dephosphorylation
  function x = P(S);
    x = k_RyR_PKA*C_cyt/(K_RyR_PKA+RyR_md*S);
  end
  function x = DP(S);
    x = k_RyR_PP*PP_cyt/(K_RyR_PP+RyR_md*S);
  end

  %Precalculate phosphorylation and dephosphorylation
  o1_o1p = f_RyR*P(o1);
  o1p_o1 = f_RyR*k_on_a*k_off_ap/(k_on_ap*k_off_a)*DP(o1p);
  o2_o2p = f_RyR*P(o2);
  o2p_o2 = f_RyR*k_on_a*k_off_ap*k_on_b*k_off_bp/(k_on_ap*k_off_a*k_on_bp*k_off_b)*DP(o2p);
  c1_c1p = P(c1);
  c1p_c1 = DP(c1p);
  c2_c2p = f_RyR*P(c2);
  c2p_c2 = f_RyR*k_on_a*k_off_ap*k_on_c*k_off_cp/(k_on_ap*k_off_a*k_on_cp*k_off_c)*DP(c2p);

  %Prepare MSM Transition matrix rows
  O1 =  [-(k_on_b*Ca_md^m+k_on_c+o1_o1p+k_off_a),k_off_b,k_on_a*Ca_md^n,k_off_c,o1p_o1,0,0,0];
  O2 =  [k_on_b*Ca_md^m,-(o2_o2p+k_off_b),0,0,0,o2p_o2,0,0];
  C1 =  [k_off_a,0,-(k_on_a*Ca_md^n+c1_c1p),0,0,0,c1p_c1,0];
  C2 =  [k_on_c,0,0,-(c2_c2p+k_off_c),0,0,0,c2p_c2];
  O1p = [o1_o1p,0,0,0,-(o1p_o1+k_on_bp*Ca_md^m+k_on_cp+k_off_ap),k_off_bp,k_on_ap*Ca_md^n,k_off_cp];
  O2p = [0,o2_o2p,0,0,k_on_bp*Ca_md^m,-(o2p_o2+k_off_bp),0,0];
  C1p = [0,0,c1_c1p,0,k_off_ap,0,-(c1p_c1+k_on_ap*Ca_md^n),0];
  C2p = [0,0,0,c2_c2p,k_on_cp,0,0,-(c2p_c2+k_off_cp)];

  %Assemble Rows
  Q = [O1;O2;C1;C2;
       O1p;O2p;C1p;C2p];

  %Calculate new states
  dS_RyR_pre = Q*Y;
  dS_RyR = dS_RyR_pre([1,2,4,5,6,7,8]);
endfunction
