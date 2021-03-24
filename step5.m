function [eta,rho_es] = step5(q_mk,ln_alpha_c_es,alpha_c_es,X_2_es,...
            ln_omega_es,ln_1minus_omega,ln_alpha_ind_es)
        
global L M
eta = zeros(L,M);

for m = 1:M
    for l = 1:L
        eta(l,m) = sum(q_mk(m,:).*ln_alpha_c_es(l,:)-alpha_c_es(l,:)*X_2_es(l,m))...
            + ln_omega_es - ln_1minus_omega - ln_alpha_ind_es(l,m) + 1;
    end
end

rho_es = 1./(1+exp(-eta)) ;

end

