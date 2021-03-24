function [alpha_ind_es,ln_alpha_ind_es] = step3(rho_es,X_2_es)

global a b
a_ind_es = a + 1 - rho_es;
b_ind_es = b + (1-rho_es).*X_2_es;
alpha_ind_es = a_ind_es./b_ind_es;
ln_alpha_ind_es = psi(a_ind_es) - log(b_ind_es);

end

