function [q_mk] = step9(rho_es,ln_alpha_c_es...
    ,alpha_c_es,X_2_es,ln_pi_es,ln_1minus_pi_es)

global L M 

phi = zeros(L,M);
phi_1 = zeros(L,M);
phi_2 = zeros(L,M);
q_mk = zeros(M,M);

for m = 1:M
    for k = 1:M
        phi_1(m,k) = sum(rho_es(:,m).*ln_alpha_c_es(:,k));
        phi_2(m,k) = sum(rho_es(:,m).*alpha_c_es(:,k).*X_2_es(:,m));
        phi(m,k) = phi_1(m,k) - phi_2(m,k) + ln_pi_es(k) +sum(ln_1minus_pi_es(1:1-k));
    end
end

exp_phi = exp(phi);
for m = 1:M
    for k = 1:M
        q_mk(m,k) = exp_phi(m,k)/sum(exp_phi(m,:));
    end
end

end

