clear all
syms f mult m w real positive
syms r_k k r_b b real

% define helper ratios r_b and r_k
% r_b = (b - B_p)/B_p = B_h/B_p and r_k = (k - K_p)/K_p = K_h/K_p
B_p = b ./ (1 + r_b); 
K_p = k ./ (1 + r_k); 
B_h = r_b .* B_p;
K_h = r_k .* K_p;

% transfer function definitions
b_sat = B_h + mult.*B_p;
k_sat = K_h + mult.*K_p;
F_over_X_sat   = sqrt( (mult.*K_p).^2 + (mult.*B_p.*w).^2);
F_over_X_unsat = sqrt( (K_p).^2       + (B_p.*w).^2);
X_sat        = 1/sqrt( (b_sat.*w).^2       + (k_sat-m.*w.^2).^2);
X_unsat      = 1/sqrt( (b.*w).^2           + (k-m.*w.^2).^2);

% solve for saturated and unsaturated force
F_sat = F_over_X_sat .* X_sat;
F_unsat = F_over_X_unsat .* X_unsat;
F_ratio = F_sat ./ F_unsat;

% find coefficents to express F_ratio in terms of mult
% form: f^2 = (cn .* [mult2,mult,1]) / (cd .* [mult2,mult,1])
[num,den] = numden(F_ratio^2);
cn = coeffs(num,mult,'all'); % cn: coeffs of numerator
cd = coeffs(den,mult,'all'); % cd: coeffs of denominator

% rearrange above form to get (cn - f^2 * cd) * [mult2,mult,1] = 0
quad = cn - f^2 * cd;
a_q = quad(1)
b_q = quad(2)
c_q = quad(3)

filename = 'simulation/modules/symbolic/get_abc_symbolic';
matlabFunction(a_q,b_q,c_q,'File',filename,'Vars',[f,m,b,k,w,r_b,r_k])
