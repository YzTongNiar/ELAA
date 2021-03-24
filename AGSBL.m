%% Initialization

% global M N L Np a b c e d f g h It Et
% global lamda omega beta rho alpha_ind alpha_c pai w
global M L N Np a b c d e f g h A Y It Et

M = 100; % Number of antenna
L = 144; % Number of multipaths
N = 2048; % Number of subcarriers
Np = 144; % Number of pilots

a = 1e-2; % Predefine parameters
b = 1e-6;
c = 2e-2;
d = 1e-2;
e = 1e-4;
f = 1e-6;
g = 0.8;
h = 0.2;
It = 30; % Maximum iterations
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
tau_1 = zeros(100,1);
for k = 1:M-1
    tau_1(k) = sum(q_mk(:,k))+1;
end

tau_2 = zeros(100,1);
for i = 1:M-1
    tau_2(i) = sum(sum(q_mk(:,i+1:M)))+c/d;
end

q_pi = zeros(M,1);
for i=1:M-1
    q_pi(i) = betarnd(tau_1(i),tau_2(i));
end
q_pi(M) = 1; 

% Channel Matirx
[interval,alpha_dp] = DP(w,alpha_c,M,L);
X = channel(alpha_ind,alpha_dp,rho,M,L);

% Obervation
% DFT matrix Np*L
A = dftmtx(N) ;
A = A(1:Np,1:L);

% generate pilots
theta = zeros(Np,1);
for i = 1:Np
    theta(i) = unifrnd(0,2*pi);
end
for i = 1:Np
    theta(i) = exp(1i*theta(i));
end 
p = diag(theta);

% matreix A Np*L
A = p*A;

% generate complex Gaussian noise Np*M
n = zeros(Np,M);
for i = 1:Np
    for j = 1:M
        n(i,j) = normrnd(0,sqrt(1/beta))+1i*normrnd(0,sqrt(1/beta));
    end
end
n = n/sqrt(2);

% obserations Y Np*M
Y = A*X+n;

%% Estimation Initialization
% set initial channel estimation as all zeros
X_es = zeros(L,M); % store the estimation of the last iteration
mu = zeros(L,M);   % store the estiamtion of the current iteration 
X_2_es = zeros(L,M); % store the estimatuion of norm |X|^2

% set initial parameter estimation
beta_es = e/f;
omega_es = g/(g+h);
lamda_es = c/d;
ln_omega_es = psi(g) - psi(g+h);
ln_1minus_omega_es = psi(h) - psi(g+h);
rho_es = g/(g+h)*ones(L,M); 
alpha_c_es = a/b*ones(L,M);
alpha_ind_es = a/b*ones(L,M);

% set iteration index
it = 0;

%% While Loop
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
    disp(it +" : "+norm(X_es-X,2));
end
