function pulse = pulse(t, t_start, dur, amp, k)
    %Produces a smoothed out stimulation current pulse.
    t_rel = t - t_start;

    rise = 1 / (1 + e^(-k*t_rel));
    fall = 1 / (1 + e^(k*(t_rel-dur)));

    pulse = amp*rise*fall;
endfunction
