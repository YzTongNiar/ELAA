function [ln_pi_es,ln_1minus_pi_es] = step7(q_mk,lamda_es)

global M
tau_1 = zeros(M,1);
tau_2 = zeros(M,1);

for k = 1:M
tau_1(k) = sum(q_mk(:,k)) + 1;   
end

for k = 1:M
tau_2(k) = sum(sum(q_mk(:,k+1:M))) + lamda_es;
end

ln_pi_es = psi(tau_1) - psi(tau_1+tau_2);
ln_1minus_pi_es = psi(tau_2)-psi(tau_1+tau_2);

end