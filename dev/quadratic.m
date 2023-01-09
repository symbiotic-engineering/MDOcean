clear all
syms f mult m w real positive
syms r_k k r_b b real

% define helper ratios r_b and r_k
% r_b = b/B_p = (B_h+B_p)/B_p and r_k = k/K_p = (K_h+K_p)/K_p
B_p = b ./ r_b; 
K_p = k ./ r_k; 
B_h = (r_b-1) .* B_p;
K_h = (r_k-1) .* K_p;

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

% find coefficents to express F_ratio (f) in terms of mult
% form: f^2 = (cn .* [mult^2,mult,1]) / (cd .* [mult^2,mult,1])
[num,den] = numden(F_ratio^2);
cn = coeffs(num,mult,'all'); % cn: coeffs of numerator
cd = coeffs(den,mult,'all'); % cd: coeffs of denominator

% rearrange above form to get (cn - f^2 * cd) * [mult2,mult,1] = 0
quad = cn - f^2 * cd;
a_q = quad(1)
b_q = quad(2)
c_q = quad(3)

bases = [b^2*w^2*f^2, ...
         b^2*w^2*f^2 , ...
         (k^2 + m^2*w^4 + b^2*w^2 - 2*k*m*w^2)*f^2, ...
         (k*m*w^2 - k^2)*f^2 , ...
         k^2*f^2 ];

for abc=1:length(quad)
matrix_tmp = coeffs(quad(abc),[r_k r_b],'all');
matrix(abc,:) = matrix_tmp([7 4 1 2 3]);
end
matrix = matrix ./ bases;
matrix = simplify(matrix,'IgnoreAnalyticConstraints',true)

latex(matrix)
latex(bases .* [r_k^2, r_k^2*r_b, r_k^2*r_b^2, r_k*r_b^2, r_b^2])

% verify that setting m=1 (so a*1^2 + b*1 + c = a+b+c = 0) causes f=1
eqn_m_1 = sum(quad) == 0;
soln_f = solve(eqn_m_1,f);
assert(soln_f == 1)

filename = 'mdocean/simulation/modules/dynamics/get_abc_symbolic';
matlabFunction(a_q,b_q,c_q,'File',filename,'Vars',[f,m,b,k,w,r_b,r_k]);
