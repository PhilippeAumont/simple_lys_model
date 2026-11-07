function dS_Na = Fast_Na(S, C_cav, V)
  %Unpack state vector
  o = S(1);
  c1 = S(2);
  c2 = S(3);
  if_na = S(4);
  i1 = S(5);
  i2 = S(6);
  ic2 = S(7);
  ic3 = S(8);
  op = S(9);
  c1p = S(10);
  c2p = S(11);
  c3p = S(12);
  ifp_na = S(13);
  i1p = S(14);
  i2p = S(15);
  ic2p = S(16);
  ic3p = S(17);

  c3 = 1 - (o+c1+c2+if_na+i1+i2+ic2+ic3+op+c1p+c2p+c3p+ifp_na+i1p+i2p+ic2p+ic3p);

  Y = [o;c1;c2;c3;if_na;i1;i2;ic2;ic3;op;c1p;c2p;c3p;ifp_na;i1p;i2p;ic2p;ic3p];

  %Parameters
  k_PKA = 6.8400e-6;   %1/uM ms
  k_PP = 1.9804e-5;    %1/uM ms
  K_PKA = 5.49415e-3;  %uM
  K_PP = 0.393025;     %uM
  PP = 0.2;              %uM

  %Calculations
  a11 = 3.802/(0.1027*e^(-(V-2.5)/17.0) + 0.20*e^(-(V-2.5)/150.0));
  a12 = 3.802/(0.1027*e^(-(V-2.5)/15.0) + 0.23*e^(-(V-2.5)/150.0));
  a13 = 3.802/(0.1027*e^(-(V-2.5)/12.0) + 0.25*e^(-(V-2.5)/150.0));
  B11 = 0.1917*e^(-(V-2.5)/20.3);
  B12 = 0.20*e^(-(V-7.5)/20.3);
  B13 = 0.22*e^(-(V-12.5)/20.3);
  a3 = 7.0e-7*e^(-(V+7.0)/7.7);
  B3 = 0.0084 + 0.00002*(V+7.0);
  a2 = 1.0/(0.188495*e^(-(V+7.0)/16.6) + 0.393956);
  B2 = a13*a2*a3/(B13*B3);
  a4 = a2/100.0;
  B4 = a3;
  a5 = a2/95000;
  B5 = a3/50.0;

  %Functions to calculate phosphorylation and dephosphorylation rates
  function x = P(S);
    x = k_PKA*C_cav/(K_PKA+S);
  end
  function x = DP(S);
    x = k_PP*PP/(K_PP+S);
  end

  %Pre-calculate phospho/dephospho rates
  o_op = P(o);
  op_o = DP(op);
  c1_c1p = P(c1);
  c1p_c1 = DP(c1p);
  c2_c2p = P(c2);
  c2p_c2 = DP(c2p);
  c3_c3p = P(c3);
  c3p_c3 = DP(c3p);
  if_ifp = P(if_na);
  ifp_if = DP(ifp_na);
  i1_i1p = P(i1);
  i1p_i1 = DP(i1p);
  i2_i2p = P(i2);
  i2p_i2 = DP(i2p);
  ic2_ic2p = P(ic2);
  ic2p_ic2 = DP(ic2p);
  ic3_ic3p = P(ic3);
  ic3p_ic3 = DP(ic3p);

  %Prepare MSM Transition matrix rows
  O =     [-(o_op+B13+a2),a13,0,0,B2,0,0,0,0,op_o,0,0,0,0,0,0,0,0];
  C1 =    [B13,-(a13+c1_c1p+B12+B3),a12,0,a3,0,0,0,0,0,c1p_c1,0,0,0,0,0,0,0];
  C2 =    [0,B12,-(a12+c2_c2p+B11+B3),a11,0,0,0,a3,0,0,0,c2p_c2,0,0,0,0,0,0];
  C3 =    [0,0,B11,-(a11+c3_c3p+B3),0,0,0,0,a3,0,0,0,c3p_c3,0,0,0,0,0];
  IF_Na = [a2,B3,0,0,-(a4+B2+a3+if_ifp+B12),B4,0,a12,0,0,0,0,0,ifp_if,0,0,0,0];
  I1 =    [0,0,0,0,a4,-(a5+i1_i1p+B4),B5,0,0,0,0,0,0,0,i1p_i1,0,0,0];
  I2 =    [0,0,0,0,0,a5,-(i2_i2p+B5),0,0,0,0,0,0,0,0,i2p_i2,0,0];
  IC2 =   [0,0,B3,0,B12,0,0,-(a12+a3+ic2_ic2p+B11),a11,0,0,0,0,0,0,0,ic2p_ic2,0];
  IC3 =   [0,0,0,B3,0,0,0,B11,-(a11+a3+ic3_ic3p),0,0,0,0,0,0,0,0,ic3p_ic3];
  Op =    [o_op,0,0,0,0,0,0,0,0,-(op_o+B13+a2),a13,0,0,B2,0,0,0,0];
  C1p =   [0,c1_c1p,0,0,0,0,0,0,0,B13,-(c1p_c1+a13+B12+B3),a12,0,a3,0,0,0,0];
  C2p =   [0,0,c2_c2p,0,0,0,0,0,0,0,B12,-(c2p_c2+a12+B11+B3),a11,0,0,0,a3,0];
  C3p =   [0,0,0,c3_c3p,0,0,0,0,0,0,0,B11,-(c3p_c3+a11+B3),0,0,0,0,a3];
  IFp_Na =[0,0,0,0,if_ifp,0,0,0,0,a2,B3,0,0,-(ifp_if+a4+B2+a3+B12),B4,0,a12,0];
  I1p =   [0,0,0,0,0,i1_i1p,0,0,0,0,0,0,0,a4,-(i1p_i1+a5+B4),B5,0,0];
  I2p =   [0,0,0,0,0,0,i2_i2p,0,0,0,0,0,0,0,a5,-(i2p_i2+B5),0,0];
  IC2p =  [0,0,0,0,0,0,0,ic2_ic2p,0,0,0,B3,0,B12,0,0,-(ic2p_ic2+a12+a3+B11),a11];
  IC3p =  [0,0,0,0,0,0,0,0,ic3_ic3p,0,0,0,B3,0,0,0,B11,-(ic3p_ic3+a11+a3)];

  %Assemble Rows
  Q = [O;C1;C2;C3;IF_Na;I1;I2;IC2;IC3;
       Op;C1p;C2p;C3p;IFp_Na;I1p;I2p;IC2p;IC3p];

  %Calculate new states
  dS_Na_pre = Q*Y;
  dS_Na = dS_Na_pre([1,2,3,5:18]);
endfunction
