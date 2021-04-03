function [group_rate] = CheckGroup(alpha_es)
% check group successful rate
% global M
% common = zeros(15,M);
% for i = 1:4
%     common(:,1+(i-1)*25:25+(i-1)*25) = alpha_es(41+(i-1)*15:55+(i-1)*15,1+(i-1)*25:25+(i-1)*25);
% end
% non_sparse = 0;
% for m = 1:M
%     for l = 1:15
%         if (common(l,m)) > 100
%             non_sparse = non_sparse + 1;
%             break
%         end
%     end
% end
% group_rate = non_sparse/100;
common = zeros(2,100);
for j = 1:4
common(1,1+(j-1)*25:25+(j-1)*25) = alpha_es(40+(j-1)*15,1+(j-1)*25:25+(j-1)*25);
common(2,1+(j-1)*25:25+(j-1)*25) = alpha_es(56+(j-1)*15,1+(j-1)*25:25+(j-1)*25);
end

count = 0;
for i = 1:100
    if (common(1,i)==0 && common(2,i)==0)
        count = count+1;
    end
end
group_rate = 1-count/100;    
end

