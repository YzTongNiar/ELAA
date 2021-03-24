function [ln_omega_es,ln_1minus_omega] = step6(rho_es)

global g h

g_es = g + sum(sum(rho_es));
h_es = h + sum(sum(1-rho_es));
ln_omega_es = psi(g_es)+psi(g_es+h_es);
ln_1minus_omega = psi(h_es)-psi(g_es+h_es);

end

