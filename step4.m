function [betta_es] = step4(mu,T)

global M Np e f Y A
e_es = M*Np+e;
f_es = f;

for i = 1:M
    f_es = f_es + norm(Y(:,i)-A*mu(:,i),2)+trace(A*T{i}*A');
end

betta_es = e_es/f_es;

end

