function set_starting_conditions()
  %Set up lysosome model
  init_Aeff = 0.30;
  init_Ca_T = 6e-3*1e6; %uM
  init_Cl = 0.001*1e6; %uM
  init_H = 0; %, proton differential from start
  init_K = 0.05*1e6; %uM
  init_Na = 0.02*1e6; %uM
  init_pH = 5;

  X0_lys = [init_H; init_pH; init_Ca_T];

  %Initial Markov State conditions
  S_LCC_cav_0 = [
  0.320206e-11,   %O
  0.524483e-2,    %C2
  0.105944e-4,    %C3
  0.951124e-8,    %C4
  0.320207e-11,   %Cp
  0.308577e-11,   %I1
  0.217536e-7,    %I2
  0.209641e-7,    %I3
  0.562222e-10,   %O-p
  0.206347e-1,    %C1-p
  0.421668e-3,    %C2-p
  0.323128e-5,    %C3-p
  0.110051e-7,    %C4-p
  0.140555e-10,   %Cp-p
  0.541817e-10,   %I1-p
  0.100683e-6,    %I2-p
  0.970287e-7     %I3-p
  ];

  S_LCC_ecav_0 = [
  0.286851e-11,   %O
  0.469850e-2,    %C2
  0.949082e-5,    %C3
  0.852050e-8,    %C4
  0.286852e-11,   %Cp
  0.276420e-11,   %I1
  0.194870e-7,    %I2
  0.187798e-7     %I3
  0.328449e-9,    %O-p
  0.120548,        %C1-p
  0.246338e-2,    %C2-p
  0.188771e-4,    %C3-p
  0.642918e-7,    %C4-p
  0.821123e-10,   %Cp-p
  0.316528e-9,    %I1-p
  0.588189e-6,    %I2-p
  0.566840e-6     %I3-p
  ];

  S_RyR_0 = [
  0.854737e-5,    %O1
  0.360412e-10,   %O2
  0.961561e-4,    %C2
  0.526065e-7,    %O1-p
  0.369705e-12,   %O2-p
  0.367832e-2,    %C1-p
  0.986431e-6     %C2-p
  ];

  S_Na_0 = [
  0.367777e-6,    %O
  0.161178e-3,    %C1
  0.132248e-1,    %C2
  0.153271e-3,    %IF
  0.146044e-4,    %I1
  0.545874e-7,    %I2
  0.125760e-1,    %IC2
  0.414822,        %IC3
  0.515006e-7,    %O-p
  0.225696e-4,    %C1-p
  0.185179e-2,    %C2-p
  0.610809e-1,    %C3-p
  0.214630e-4,    %IF-p
  0.217162e-5,    %I1-p
  0.301835e-7,    %I2-p
  0.176099e-2,    %IC2-p
  0.580859e-1     %IC3-p
  ];

  S_IKr_0 = [
  0.332600e-3,    %O
  0.135218e-2,    %C1
  0.873596e-3,    %C2
  0.763767e-4     %I
  ];

  S_RyR_md_0 = [
  0.854737e-5,    %O1
  0.360412e-10,   %O2
  0.961561e-4,    %C2
  0.526065e-7,    %O1-p
  0.369705e-12,   %O2-p
  0.367832e-2,    %C1-p
  0.986431e-6     %C2-p
  ];

  %Initial Conditions
  ans = [
  -78.2787,       %membrane potential
  0.100157,       %myoplasmic Ca (Ca_i)
  0.100157,       %subspace Ca (Ca_ss)
  0.100157,       %microdomain Ca (Ca_md)
  1081.23,        %JSR Ca (Ca_JSR)
  1081.23,        %NSR Ca (Ca_NSR)
  8.66981,        %LTRPNCa
  123.369,        %HTRPNCa
  10508.5,        %Myoplasmic Na
  145400,         %Myoplasmic K
  0.533799e-2,   %a_to_f
  0.999945,       %i_to_f
  0.713943e-3,   %a_ur
  0.996991,       %i_ur
  0.713943e-3,   %a_Kss
  0.225905,       %f_cav_PLM_p
  0.908852,       %f_ecav_IKur
  0.713943e-3,   %a_urp
  0.996991,       %i_urp
  0.252661,       %f_ecav_IKto,f
  0.111499e-2,   %a_to_fp
  0.999983,       %i_to_fp
  0.186637,       %f_cyt_PLB_p
  0.364102,       %f_cyt_Tnl_p
  0.799452e-3,   %R_cav_PKA
  0.626341e-27,  %R_cav_GRK2
  0.132189e-2,    %Gs_cav_aGTP
  0.180824e-2,    %Gs_cav_By
  0.487356e-3,    %Gs_cav_aGDP
  0.478002e-1,    %R_ecav_PKA
  0.626341e-27,   %R_ecav_GRK2
  0.230801e-1,    %Gs_ecav_aGTP
  0.237276e-1,    %Gs_ecav_By
  0.648475e-3,    %Gs_ecav_aGDP
  0.155949e-2,    %R_cyt_PKA
  0.626341e-27,   %R_cyt_GRK2
  0.331511e-3,    %Gs_cyt_aGTP
  0.663570e-3,    %Gs_cyt_By
  0.333058e-3,    %Gs_cyt_aGDP
  0.000000,        %cAMP_cav_AC56
  0.000000,        %cAMP_ecav_AC47
  0.000000,        %cAMP_cyt_AC56
  0.000000,        %cAMP_cyt_AC47
  0.125103e-1,    %PDE3_cav_p
  0.580798e-2,    %PDE4_cav_p
  0.000000,        %cAMP_cav_PDE2
  0.000000,        %cAMP_cav_PDE3
  0.000000,        %cAMP_cav_PDE4
  0.158226e-1,    %PDE4_ecav_p
  0.000000,        %cAMP_ecav_PDE2
  0.000000,        %cAMP_ecav_PDE4
  0.120998e-2,    %PDE3_cyt_p
  0.373102e-2,    %PDE4_cyt_p
  0.000000,        %cAMP_cyt_PDE2
  0.000000,        %cAMP_cyt_PDE3
  0.000000,        %cAMP_cyt_PDE4
  7.92317,         %cAMP_cav_PKA
  0.299288,        %ARC_cav
  0.303358e-1,    %A2RC_cav
  0.858440,        %A2R_cav
  0.459397e-1,    %C_cav
  0.823499,        %PKIC_cav
  6.74029,         %cAMP_ecav_PKA
  0.653988,        %ARC_ecav
  0.132861,        %A2RC_ecav
  1.17000,         %A2R_ecav
  0.147623,        %C_ecav
  1.03338,         %PKIC_ecav
  9.32461,         %cAMP_cyt_PKA
  0.996350e-1,    %ARC_cyt
  0.140099e-1,    %A2RC_cyt
  0.273868,        %A2R_cyt
  0.665022e-1,    %C_cyt
  0.218365,        %PKIC_cyt
  0.213571e-1,    %Inhib1_cyt_p
  0.253399,        %cAMP_cav
  0.507889,        %cAMP_ecav
  0.407775,        %cAMP_cyt
  0.254152e-11,   %RyR Modulation factor
  S_LCC_cav_0,
  S_LCC_ecav_0,
  S_RyR_0,
  S_RyR_md_0,
  S_Na_0,
  S_IKr_0,
  X0_lys
  ];

  save("X0_generated.mat", "ans");
endfunction
