% corresponding to 10,15,20,25,30 dB
beta_list = [10 10^1.5 100 10^2.5 1000];
err_list = zeros(1,5);
err_list2 = zeros(1,5);

% LS & genie LS estimation
for i = 1:size(beta_list,2)
    EsInit;
    n = CNormal(0,1/beta_list(i),L,M);
    Y = A*X + n;
    X_es_lse = A\Y;
    X_es_gls = DeleteSparse(A\Y);
    err_list(i) = norm(X_es_lse-X,2)/(L*M);
    err_list2(i) = norm(X_es_gls-X,2)/(L*M);
end

x = 10:5:30;
%plot(x,err_list,x,err_list2,x,err_list3);
plot(x,err_list,x,err_list2,x,err_list3);
set(gca,'xtick',10:5:30);
xlabel('SNR/dB');
ylabel('MSE');
legend('LS','genie LS','AGSBL')
hold on;