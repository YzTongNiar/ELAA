
% step1 get estimation of q(x)
d_vector = g/(g+h)*a/b*sum(q_mk,2)+(1-g/(g+h)*a/b);
d_matrix = ones(L,M);
for i = 1:M
    d_matrix(:,i) = d_matrix(:,i)*d_vector(i);
end

mu = zeros(L,M);
T = cell(100,1);
for i = 1:M
    T{i} = inv(e/f*(A')*A+diag(d_matrix(:,i)));
    mu(:,i) = e/f*T{i}*(A')*Y(:,i);
end

% step2 get estimation of q(alpha_c)

a_kl = zeros(L,M);
row_sum = sum(q_mk,1);
for k = 1:M
    for l = 1:L
        a_kl(l,k)=a+g/(g+h)*row_sum(k);
    end
end
        

X_2 = zeros(L,M);
for i = 1:L
    for j = 1:M
        T_m = T{j};
        X_2(i,j) = mu(i,j)^2+T_m(i,i);
    end
end
b_kl = (a_kl-a).*X_2 + b;
alpha_c_es = a_kl./b_kl;
In_alpha_c_es = psi(a_kl) - log(b_kl);


 a_kl = a+rho_es(1,:)*q_mk(:,2);
 a_kl_2 = a + sum(rho_es(1,:).*q_mk(:,2)');
