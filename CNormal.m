function [X] = CNormal(mu,variance,L,M)
% generate a complex gaussian matrix
X = normrnd(mu,sqrt(variance/2),[L M]) + 1i*normrnd(0,sqrt(variance/2),[L M]);
end

