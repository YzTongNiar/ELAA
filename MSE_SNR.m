beta_list = [10 10^1.5 100 10^2.5 1000];
err_list = zeros(3,5);
for i = 1:size(beta_list,2)
    beta_es = beta_list(i);
    n = CNormal(0,1/beta_list(i),144,100);
    Y = A*X+n;
    it = 0;
    % First Do LS & genie-LS Estimation
    X_es_lse = A\Y;
    X_es_gls = DeleteSparse(A\Y);
    err_list(1,i) = norm(X_es_lse-X,2)/(L*M);
    err_list(2,i) = norm(X_es_gls-X,2)/(L*M);
while(1)
    % step1
    [mu,T] = step1(beta_es,rho_es,alpha_c_es,alpha_ind_es,q_mk);
    % step2
    [alpha_c_es,ln_alpha_c_es,X_2_es] = step2(rho_es,mu,q_mk,T);
    % step3
    [alpha_ind_es,ln_alpha_ind_es] = step3(rho_es,X_2_es);
    % step4
    [beta_es] = step4(mu,T);
    % step5
    [eta,rho_es] = step5(q_mk,ln_alpha_c_es,alpha_c_es,X_2_es,...
            omega_es,ln_1minus_omega_es,ln_alpha_ind_es);
    % step6
    [ln_omega_es,ln_1minus_omega_es] = step6(rho_es);
    % step7
    [ln_pi_es,ln_1minus_pi_es] = step7(q_mk,lamda_es);
    % step8
    [lamda_es] = step8(ln_1minus_pi_es);
    % step9
    [q_mk] = step9(rho_es,ln_alpha_c_es...
    ,alpha_c_es,X_2_es,ln_pi_es,ln_1minus_pi_es);
    % step10
    it = it+1;    
    % while loop condition
    if (it>It||check_rel(X_es,mu)) % X_es is from previos iteration; m is from current iteration
        break;
    end    
    % Update the estiamtion
    X_es = mu;
end
%err_list3(i) = norm(X_es-X,2)/(L*M);
err_list(3,i) = norm(X_es-X,2)/(L*M);
end
x = 10:5:30;
%plot(x,err_list,x,err_list2,x,err_list3);
plot(x,err_list(1,:),x,err_list(2,:),x,err_list(3,:));
set(gca,'xtick',10:5:30);
xlabel('SNR/dB');
ylabel('MSE');
legend('LS','genie LS','AGSBL')
hold on;