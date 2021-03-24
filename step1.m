function [mu,T] = step1(betta_es,rho_es,alpha_c_es,alpha_ind_es,q_mk)

global M L A Y
temp = zeros(L,M);
d_matrix = zeros(L,M);

for l = 1:L
    for m =1:M
        temp(l,m) = sum(q_mk(m,:).*alpha_c_es(l,:));
    end
end
        
for l = 1:L
    for m =1:M
        d_matrix(l,m) = (1-rho_es(l,m))*alpha_ind_es(l,m)...
            + rho_es(l,m)*temp(l,m);
    end
end

mu = zeros(L,M);
T = cell(M,1);

for i = 1:M
    T{i} = inv(betta_es*(A')*A+diag(d_matrix(:,i)));
    mu(:,i) = betta_es*T{i}*(A')*Y(:,i);
end

end

