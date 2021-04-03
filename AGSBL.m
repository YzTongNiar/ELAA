%% Initialization

global M L N Np a b c d e f g h A X Y n It Et beta IndIndex

M = 100; % Number of antenna
L = 144; % Number of multipaths
N = 2040; % Number of subcarriers; Not implemented yet
Np = 144; % Number of pilots; Not implemented yet

a = 1e-2; % Predefine parameters
b = 1e-6;
c = 2e-2;
d = 1e-2;
e = 1e-4;
f = 1e-6;
g = 0.8;
h = 0.2;
It = 100; % Maximum iterations
Et = 1e-4; % Error tolerance

% random parameter generating
% variance of noise
beta = gamrnd(1/f,e);

% generate hybird parameter
omega = betarnd(g,h);
rho = zeros(L,M);

for i=1:L
    for j = 1:M
        rho(i,j) = binornd(1,omega);
    end
end

% generae the individual and common sparsity parameters
alpha_ind = gamrnd(1/b,a,L,M); % correspond to alpha_c_m,l
alpha_c = gamrnd(1/b,a,L,M); % correspond to alpha*_k,l

% generate DP parameters
lamda = gamrnd(1/d,c);

pai = betarnd(1,lamda,M,1);
pai(M) = 1;
% the probability for each atoms using stick breaking
w = zeros(M,1);
for i = 1:M
    w(i) = pai(i);
    w_2 = 1;
    for j = 1:i-1
        w_2 = w_2*(1-pai(j));
    end
    w(i) = w(i)*w_2;
end

% the probability alpha_c(:,m) is from alpha_dp(:,k)
q_mk = zeros(M,M);
for i = 1:M
    q_mk(i,:) = w';
end

% using quation(36) to get q_pi
tau_1 = zeros(M,1);
for k = 1:M-1
    tau_1(k) = sum(q_mk(:,k))+1;
end

tau_2 = zeros(M,1);
for i = 1:M-1
    tau_2(i) = sum(sum(q_mk(:,i+1:M)))+c/d;
end

q_pi = zeros(M,1);
for i=1:M-1
    q_pi(i) = betarnd(tau_1(i),tau_2(i));
end
q_pi(M) = 1; 

% Channel Matirx
[X, IndIndex] = channel2(L,M);

% Obervation
%DFT matrix Np*L
A = dftmtx(L) ; 

% generate noise
n = CNormal(0,1/beta,144,100);

% obserations Y Np*M
Y = A*X+n;

%% Estimation Initialization
% set initial channel estimation as all zeros
X_es = zeros(L,M); % store the estimation of the last iteration
mu = zeros(L,M);   % store the estiamtion of the current iteration 
X_2_es = zeros(L,M); % store the estimatuion of norm |X|^2

% set initial parameter estimation
%beta_es = e/f;
beta_es = 10;
omega_es = g/(g+h);
lamda_es = c/d;
ln_omega_es = psi(g) - psi(g+h);
ln_1minus_omega_es = psi(h) - psi(g+h);
rho_es = g/(g+h)*ones(L,M); 
alpha_c_es = a/b*ones(L,M);
alpha_ind_es = a/b*ones(L,M);

% set arrays to store the data
err_list3 = zeros(1,5);
beta_list = [10 10^1.5 100 10^2.5 1000];

%% While Loop
for i = 1:size(beta_list,2)
    beta = beta_list(i);
    beta_es = beta_list(i);
    n = CNormal(0,1/beta,144,100);
    Y = A*X+n;
    it = 0;
while(1)
    % step1
    [mu,T] = step1(beta_es,rho_es,alpha_c_es,alpha_ind_es,q_mk);
    disp("step1 finish");
    % step2
    [alpha_c_es,ln_alpha_c_es,X_2_es] = step2(rho_es,mu,q_mk,T);
     disp("step2 finish");
    % step3
    [alpha_ind_es,ln_alpha_ind_es] = step3(rho_es,X_2_es);
     disp("step3 finish");
    % step4
    [beta_es] = step4(mu,T);
     disp("step4 finish");
    % step5
    [eta,rho_es] = step5(q_mk,ln_alpha_c_es,alpha_c_es,X_2_es,...
            omega_es,ln_1minus_omega_es,ln_alpha_ind_es);
     disp("step5 finish");
    % step6
    [ln_omega_es,ln_1minus_omega_es] = step6(rho_es);
     disp("step6 finish");
    % step7
    [ln_pi_es,ln_1minus_pi_es] = step7(q_mk,lamda_es);
     disp("step7 finish");
    % step8
    [lamda_es] = step8(ln_1minus_pi_es);
     disp("step8 finish");
    % step9
    [q_mk] = step9(rho_es,ln_alpha_c_es...
    ,alpha_c_es,X_2_es,ln_pi_es,ln_1minus_pi_es);
     disp("step9 finish");
    % step10
    it = it+1;
    
    % while loop condition
    if (it>It||check_rel(X_es,mu)) % X_es is from previos iteration; m is from current iteration
        break;
    end
    
    % Update the estiamtion
    X_es = mu;
    disp(it +" : "+norm(X_es-X,2)/(L*M));
end
err_list3(i) = norm(X_es-X,2)/(L*M);
end