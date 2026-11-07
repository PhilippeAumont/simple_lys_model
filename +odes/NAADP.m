function NAADP = NAADP(t, p)
  if strcmp(p.NAADP_protocol, "cte")
    NAADP = p.NAADP_C;

  elseif strcmp(p.NAADP_protocol, "sig")
    NAADP = p.NAADP_C/(1+exp(-p.NAADP_k1*(t-p.NAADP_t0)));

  elseif strcmp(p.NAADP_protocol, "spike")
    NAADP = odes.pulse(t, p.NAADP_stim_start, p.NAADP_stim_dur, p.NAADP_C, p.NAADP_stim_k2);
 end
endfunction
