clear all;
syms f mult m w k r_b b real positive;
syms r_k real;

% define helper ratios r_b and r_k
% r_b = b/B_p = (B_h+B_p)/B_p and r_k = k/K_p = (K_h+K_p)/K_p
B_p = b ./ r_b;
K_p = k ./ r_k;
B_h = (r_b - 1) .* B_p;
K_h = (r_k - 1) .* K_p;

% transfer function definitions
b_sat = B_h + mult .* B_p;
k_sat = K_h + mult .* K_p;
F_over_X_sat   = sqrt((mult .* K_p).^2 + (mult .* B_p .* w).^2);
F_over_X_unsat = sqrt((K_p).^2       + (B_p .* w).^2);
X_sat        = 1 / sqrt((b_sat .* w).^2       + (k_sat - m .* w.^2).^2);
X_unsat      = 1 / sqrt((b .* w).^2           + (k - m .* w.^2).^2);

% solve for saturated and unsaturated force
F_sat = F_over_X_sat .* X_sat;
F_unsat = F_over_X_unsat .* X_unsat;
F_ratio = F_sat ./ F_unsat;

% find coefficents to express F_ratio (f) in terms of mult
% form: f^2 = (cn .* [mult^2,mult,1]) / (cd .* [mult^2,mult,1])
[num, den] = numden(F_ratio^2);
cn = coeffs(num, mult, 'all'); % cn: coeffs of numerator
cd = coeffs(den, mult, 'all'); % cd: coeffs of denominator

% rearrange above form to get (cn/f^2 - cd) * [mult2,mult,1] = 0
quad = cn / f^2 - cd;
a_q = quad(1);
b_q = quad(2);
c_q = quad(3);

filename = 'mdocean/simulation/modules/dynamics/get_abc_symbolic';
matlabFunction(a_q, b_q, c_q, 'File', filename, 'Vars', [f, m, b, k, w, r_b, r_k]);

%% post processing into matrix form

% normalize to use nondimensional w/wn and zeta
syms wn zeta w_ratio real positive;
% wn^2 = k / m so m^2 = k^2/wn^4
quad = subs(quad, m^2, k^2 / wn^4);

% zeta^2 = b^2 / (4 k m) so km = b^2 / (4 zeta^2)
quad = subs(quad, k * m, b^2 / (4 * zeta^2));

% b/m = 2*zeta*wn and m=k/wn^2 so b=2 zeta k / wn
quad = subs(quad, b, 2 * zeta * k / wn);

quad = subs(quad, w / wn, w_ratio);

quad = simplify(quad / k^2); % normalize

bases = [zeta^2 * w_ratio^2, ...
         zeta^2 * w_ratio^2, ...
         (1 + w_ratio^4 + 4 * zeta^2 * w_ratio^2 - 2 * w_ratio^2), ...
         (w_ratio^2 - 1), ...
         1];

for abc = 1:length(quad)
    matrix_tmp = coeffs(quad(abc), [r_k r_b], 'all');
    matrix(abc, :) = matrix_tmp([7 4 1 2 3]);
end
matrix = matrix ./ bases;
matrix = simplify(matrix, 'IgnoreAnalyticConstraints', true);

latex(matrix);
rkb = [r_k^2, r_k^2 * r_b, r_k^2 * r_b^2, r_k * r_b^2, r_b^2];
bases_rk = bases .* rkb;
pretty(bases_rk');
latex(bases_rk);

% verify that setting m=1 (so a*1^2 + b*1 + c = a+b+c = 0) causes f=1
eqn_m_1 = sum(quad) == 0;
soln_f = solve(eqn_m_1, f);
assert(soln_f == 1);

%% replace rk with other variables
% notebook p97 2-1-25
syms omega_u_star zeta_u star real positive;

rk_new = star^2 / (star^2 - 1);
rb_new = 1 / (1 - zeta_u / (zeta * star));
w_rat_new = omega_u_star / star;

rkb_new = subs(rkb, [r_k, r_b], [rk_new, rb_new]);
pretty(rkb_new.');
bases_new = subs(bases, w_ratio, w_rat_new);
bases_rk_new = bases_new .* rkb_new;
scale = (star^2 - 1)^2 * (star * zeta)^-2 * (zeta_u - star * zeta)^2;
pretty(simplify(bases_rk_new.' * scale, 'IgnoreAnalyticConstraints', true));

%% check that parameterization used in paper equals "new" just above
% notebook p99 2-1-25
matrix_check = [-4  4 / f^2 1 / f^2 -1  0  0
                8  0     0      2  8  2
                -4 -4    -1     -1 -8 -2];
bases_unsquared_check = [omega_u_star * (zeta_u - star * zeta), ...
                         omega_u_star * star * zeta, ...
                         omega_u_star^2 - star^2, ...
                         star^2 - 1];
bases_squared_check = [bases_unsquared_check(1)^2; ...
                       bases_unsquared_check(2)^2; ...
                       bases_unsquared_check(3)^2; ...
                       bases_unsquared_check(4)^2; ...
                       bases_unsquared_check(1) * bases_unsquared_check(2); ...
                       bases_unsquared_check(3) * bases_unsquared_check(4)];

bases_should_match = [bases_unsquared_check(1)^2; ... % a^2
                      -bases_unsquared_check(1) * bases_unsquared_check(2); ... % -ab
                      bases_unsquared_check(3)^2 + 4 * bases_unsquared_check(2)^2; ... % c^2+4b^2
                      bases_unsquared_check(3) * bases_unsquared_check(4); ... % cd
                      bases_unsquared_check(4)^2]; % d^2

should_be_1 = simplify(bases_should_match ./ (bases_rk_new.' * scale));
assert(all(should_be_1 == 1));

%% and check that it equals the original (r_k, r_b, quad)
quad_new = subs(quad, [r_k, r_b, w_ratio], [rk_new, rb_new, w_rat_new]);
should_be_1_also = simplify(quad_new.' ./ (matrix_check * bases_squared_check / scale));
assert(all(should_be_1_also == 1));

%% optimal powertrain case
bases_rk = subs(bases_rk, {w_ratio, r_b}, {1, 2}).';

abc = matrix * bases_rk / 4;
pretty(abc);

discrim = sqrt(abc(2)^2 - 4 * abc(1) * abc(3));
m_plus  = simplify((-abc(2) + discrim) / (2 * abc(1)));
m_minus = simplify((-abc(2) - discrim) / (2 * abc(1)));

syms r_k_zeta;
m_plus  = subs(m_plus, r_k, r_k_zeta / zeta);
m_minus = subs(m_minus, r_k, r_k_zeta / zeta);

pretty(m_plus);
pretty(m_minus);

latex(m_plus);

m_plus_rkz_inf = simplify(limit(m_plus, r_k_zeta, Inf));

%% sub in numerical
% relationship between m, f_sat, and w for the ideal case (max power control)
% see notebook p115-118
close all;

f_sat = linspace(0.01, 1);
r_k_z = linspace(0.01, 3);

[F, RKZ] = meshgrid(f_sat, r_k_z);

mult_num = vpa(subs(m_plus, {f, r_k_zeta}, {F, RKZ}));
mult_num(:, :, 2) = vpa(subs(m_minus, {f, r_k_zeta}, {F, RKZ}));

energy = F.^2 ./ mult_num(:, :, 1);

figure;
contourf(F, RKZ, mult_num(:, :, 1));
xlabel('f_{sat}');
ylabel('|r_k \zeta|');
title('m_+');
colorbar;

figure;
contourf(F, RKZ, mult_num(:, :, 2), [min(mult_num(:, :, 2), [], 'all') 0:.1:1 max(mult_num(:, :, 2), [], 'all')]);
xlabel('f_{sat}');
ylabel('|r_k \zeta|');
title('m_-');
colorbar;
caxis([0 1]);

figure;
contourf(F, RKZ, energy);
xlabel('f_{sat}');
ylabel('|r_k \zeta|');
title('Energy scale for m_-');
colorbar;
caxis([0 1]);

r_k_z_selected = r_k_z([1 4 7 11 17 34]); % corresponds to values [.01 .1 .2 .3 .5 1]
figure;
hold on;

colors = {'m', 'r', 'y', 'g', 'c', 'b', 'k'};

for i = 1:length(r_k_z_selected)
    rkzi = r_k_z_selected(i);
    idx = RKZ == rkzi;
    f_i = F(idx);
    m_tmp = mult_num(:, :, 1);
    m_i = m_tmp(idx);
    e_i = energy(idx);

    %     subplot 131
    %     plot(f_i,m_i,colors{i})
    %     hold on
    %     xlabel('Force Saturation Factor f_{sat}')
    %     ylabel('Powertrain Control Multiplier m_{sat}')

    if i == 1 % number of decimals to round r_k_zeta in legend
        round_num = 2;
    else
        round_num = 1;
    end
    subplot 121;
    plot(f_i, e_i, colors{i});
    hold on;
    xlabel('Force Saturation Factor f_{sat}');
    ylabel('Energy Multiplier e');

    subplot 122;
    plot(m_i, e_i, colors{i}, 'DisplayName', ['|r_k \zeta| = ' num2str(round(rkzi, round_num))]);
    hold on;
    xlabel('Powertrain Control Multiplier m_{sat}');
    ylabel('Energy Multiplier e');
end

m_inf = subs(m_plus_rkz_inf, f, f_sat);
e_inf = f_sat.^2 ./ m_inf;
% subplot 131
% plot(f_sat,m_inf,'k')
subplot 121;
plot(f_sat, e_inf, 'k');
grid on;
subplot 122;
plot(m_inf, e_inf, 'k', 'DisplayName', '|r_k \zeta| \rightarrow \infty');
grid on;

leg = legend;
sgtitle('Effect of Controller Non-Idealities on Energy');
title(leg, 'Robustness Metric');
improvePlot;

%% now calculate required fsat based on r_k zeta
syms F_max_over_F_d;
% what I called F_unsat above is technically F_unsat_over_Fd
alpha = 2 / sym('pi') * (F_unsat / F_max_over_F_d * asin(F_max_over_F_d / F_unsat)) + sqrt(1 - (F_max_over_F_d / F_unsat)^2);
f_sat_new = alpha * F_max_over_F_d / F_unsat;

% do the same nondimensionalizing that I did to quad above
f_sat_new = subs(expand(f_sat_new), m^2, k^2 / wn^4);

% zeta^2 = b^2 / (4 k m) so km = b^2 / (4 zeta^2)
f_sat_new = subs(f_sat_new, k * m, b^2 / (4 * zeta^2));

% b/m = 2*zeta*wn and m=k/wn^2 so b=2 zeta k / wn
f_sat_new = subs(f_sat_new, b, 2 * zeta * k / wn);

f_sat_new = subs(f_sat_new, w / wn, w_ratio);

% ideal powertrain case
f_sat_new = subs(f_sat_new, {w_ratio, r_b}, {1, 2});
f_sat_new = subs(f_sat_new, r_k, r_k_zeta / zeta);

f_sat_new = simplify(f_sat_new);

% doesnt work - solver unable to find explicit solution. so redo numerical
% syms f_sat_var
% eqn = f_sat_var - f_sat_new == 0;
% F_max_over_F_d_soln = solve(eqn,F_max_over_F_d);
% F_max_over_F_d_num = subs(F_max_over_F_d_soln, {r_k_zeta,f_sat_var},{RKZ,F});
%
% figure
% contourf(F,F_max_over_F_d_soln,energy)

% now use F to refer to F_max/F_d instead of f_sat
m_plus_new = subs(m_plus, f, f_sat_new);

mult_new_num = vpa(subs(m_plus_new, {F_max_over_F_d, r_k_zeta, pi}, {F, RKZ, 3.14}));
FS = vpa(subs(f_sat_new, {F_max_over_F_d, r_k_zeta}, {F, RKZ}));
energy = FS.^2 ./ mult_new_num;

figure;
contourf(F, RKZ, mult_new_num);
title('Mult');
figure;
contourf(F, RKZ, energy);
title('Energy');
figure;
contourf(F, RKZ, FS);
title('f sat');
