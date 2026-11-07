function I = I_stim(t, p)
  if strcmp(p.protocol, "none")
    I = 0;

  elseif strcmp(p.protocol, "spike")
    I = odes.pulse(t, p.stim_start, p.stim_dur, p.stim_amp, p.stim_k);

  elseif strcmp(p.protocol, "train")

    t_diff = t - p.stim_start;
    t_in_cycle = mod(t_diff, p.stim_period);

    rise = 1 / (1 + exp(-p.stim_k*t_in_cycle));
    fall = 1 / (1 + exp(p.stim_k*(t_in_cycle - p.stim_dur)));

    I = p.stim_amp * rise * fall * (t_diff >= 0);

  elseif strcmp(p.protocol, "two_spikes")
    pulse1 = odes.pulse(t, p.stim_start, p.stim_dur, p.stim_amp, p.stim_k);
    pulse2 = odes.pulse(t, p.stim_2nd_start, p.stim_dur, p.stim_amp, p.stim_k);

    I = pulse1+pulse2;

 end
endfunction
