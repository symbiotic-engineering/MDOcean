clear all
syms f mult m w real positive
syms r_k k r_b b real

B_p = b ./ (1 + r_b);
K_p = k ./ (1 + r_k);
B_h = r_b .* B_p;
K_h = r_k .* K_p;
b_sat = B_h + mult.*B_p;
k_sat = K_h + mult.*K_p;
F_over_X_sat = sqrt((mult.*K_p).^2+(mult.*B_p.*w).^2);
F_over_X_unsat = sqrt((K_p).^2+(B_p.*w).^2);
X_sat = 1/sqrt((b_sat.*w).^2+(k_sat-m.*w.^2).^2);
X_unsat = 1/sqrt((b.*w).^2+(k-m.*w.^2).^2);
F_sat = F_over_X_sat .* X_sat;
F_unsat = F_over_X_unsat .* X_unsat;
F_ratio = F_sat ./ F_unsat;

[n,d] = numden(F_ratio^2);
[cn,mn] = coeffs(n,mult);
[cd,md] = coeffs(d,mult);

% cn .* mn = cd .* md * f^2;
quad = [cn 0 0] - f^2 * cd;
a_q = quad(1)
b_q = quad(2)
c_q = quad(3)

% eqn1 = collect(eqn^2,mult);
% pretty(eqn1)
% 
% syms b_q c_q real
% eqn2 = 0 == mult^2 + b_q*mult + c_q;
% soln_q = solve([eqn1,eqn2],[b_q,c_q])

%eqn = simplify(F_ratio) == f;
%soln = simplify(solve(eqn,mult))
%matlabFunction(soln)