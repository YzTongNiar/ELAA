function [alpha_c_es,In_alpha_c_es,X_2_es] = step2(rho_es,mu,q_mk,T)

global L M a b
a_kl = zeros(L,M);
b_kl = zeros(L,M);

for k = 1:M
    for l = 1:L
        a_kl(l,k)=a+rho_es(l,:)*q_mk(:,k);
    end
end
        
X_2_es = zeros(L,M);
for i = 1:L
    for j = 1:M
        T_m = T{j};
        X_2_es(i,j) = abs(mu(i,j))^2+norm(T_m(i,i));
    end
end

for l = 1:L
    for k = 1:M
        b_kl(l,k) = b + sum(rho_es(l,:).*q_mk(:,k)'.*X_2_es(l,:));
    end
end

alpha_c_es = a_kl./b_kl;
In_alpha_c_es = psi(a_kl) - log(b_kl);

end

