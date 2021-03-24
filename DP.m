function [interval,alpha_dp] = DP(w,alpha_c,M,L)

interval = ones(M,1);
alpha_dp = zeros(L,M);
for i = 1:M
    r = rand(1);
    current = w(1);
    for k = 1:M
        if r > current
            current = current + w(k+1);
            interval(i) = interval(i) + 1;
        end
    end
end
for i = 1:M
    alpha_dp(:,i) = alpha_c(:,interval(i));
end
end