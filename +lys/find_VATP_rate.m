function [Hpump_query] = find_VATP_rate(v_flux,psi_query,pH_query)
%Find proton flux for specific membrane potential and luminal pH.

pH=v_flux(1,:,2);
psi=v_flux(:,1,1);
Hpump=v_flux(:,:,3);
psi_query=real(psi_query);
pH_query = real(pH_query);
%%%for values outside the dataset's range, use the first of last point in
%%%the data set.
if psi_query > psi(end)
    psi_query = psi(end) ;
elseif psi_query < psi(1)
    psi_query = psi(1);
end

if pH_query > pH(1)
    pH_query = pH(1);
elseif pH_query < pH(end)
    pH_query = pH(end);
end
%%%for values inside the dataset's range,
Hpump_query = interp2(pH,psi,Hpump,pH_query,psi_query);

