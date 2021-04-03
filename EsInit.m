global M L a b c d g h Y
X_es = Y;
mu = zeros(L,M);   % store the estiamtion of the current iteration 
X_2_es = zeros(L,M); % store the estimatuion of |X|^2

% set initial parameter estimation
omega_es = g/(g+h);
lamda_es = c/d;
ln_omega_es = psi(g) - psi(g+h);
ln_1minus_omega_es = psi(h) - psi(g+h);
rho_es = g/(g+h)*ones(L,M); 
alpha_c_es = a/b*ones(L,M);
alpha_ind_es = a/b*ones(L,M);