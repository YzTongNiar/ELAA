function [X] = channel(alpha_ind,alpha_dp,rho,M,L)

X = zeros(L,M);

for l = 1:L
    for m = 1:M
        gauss1 = normrnd(0,sqrt(1/alpha_dp(l,m)))+1i*normrnd(0,sqrt(1/alpha_dp(l,m)));
        gauss2 = normrnd(0,sqrt(1/alpha_ind(l,m)))+1i*normrnd(0,sqrt(1/alpha_ind(l,m)));
        gauss1 = gauss1/sqrt(2);
        gauss2 = gauss2/sqrt(2);
        X(l,m) = gauss1^rho(l,m)*gauss2^(1-rho(l,m));
    end
end

end

