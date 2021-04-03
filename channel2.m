function [X,index] = channel2(L,M) % Index is the column index of individual non-zeros
X = zeros(L,M);
% wholly common non-zero entires
X(1,:) = CNormal(0,1,1,M);
% partially common non-zeros entires
for i = 1:4
    X(41+(i-1)*15:55+(i-1)*15,1+(i-1)*25:25+(i-1)*25) = CNormal(0,1,15,25);
end
% individual non-zeros entires
index = randperm(100,44);
for j = 101:144
    X(j,index(j-100)) = CNormal(0,1,1,1);
end

end

