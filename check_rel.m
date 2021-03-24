function [flag_1] = check_rel(X_es1,X_es2)

global M Et
flag_1 = 1;

for i = 1:M
    a = norm(X_es1(:,i)-X_es2(:,i),2);
    b = Et*norm(X_es1(:,i),2);
    if (norm(X_es1(:,i)-X_es2(:,i),2)> Et*norm(X_es1(:,i),2))
        flag_1 = 0;
        break;
    end
end
end

