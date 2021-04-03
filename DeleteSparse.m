function [Y_genie] = DeleteSparse(Y)
global L M IndIndex
% get rid off the sparse entires
Y_genie = zeros(L,M);
% keep the non-zero entires
Y_genie(1,:) = Y(1,:);

for i = 1:4
    Y_genie(41+(i-1)*15:55+(i-1)*15,1+(i-1)*25:25+(i-1)*25) = ...
        Y(41+(i-1)*15:55+(i-1)*15,1+(i-1)*25:25+(i-1)*25);
end

for j = 101:144
    Y_genie(j,IndIndex(j-100)) = Y(j,IndIndex(j-100));
end
end

