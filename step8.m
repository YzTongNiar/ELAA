function [lamda_es] = step8(ln_1minus_pi_es)

global M c d

c_es = c+M-1;
d_es = d - sum(ln_1minus_pi_es(1:M-1));
lamda_es = c_es/d_es;

end

