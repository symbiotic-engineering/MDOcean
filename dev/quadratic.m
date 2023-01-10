clear all
syms f mult m w k r_b b real positive
syms r_k real

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

% rearrange above form to get (cn/f^2 - cd) * [mult2,mult,1] = 0
quad = cn/f^2 - cd;
a_q = quad(1)
b_q = quad(2)
c_q = quad(3)

filename = 'mdocean/simulation/modules/dynamics/get_abc_symbolic';
matlabFunction(a_q,b_q,c_q,'File',filename,'Vars',[f,m,b,k,w,r_b,r_k]);

%% post processing into matrix form

% normalize to use nondimensional w/wn and zeta
syms wn zeta w_ratio real positive
% wn^2 = k / m so m^2 = k^2/wn^4
quad = subs(quad, m^2, k^2/wn^4);

% zeta^2 = b^2 / (4 k m) so km = b^2 / (4 zeta^2)
quad = subs(quad, k*m, b^2 / (4*zeta^2));

% b/m = 2*zeta*wn and m=k/wn^2 so b=2 zeta k / wn
quad = subs(quad, b, 2*zeta*k / wn);

quad = subs(quad, w/wn, w_ratio);

quad = simplify(quad / k^2); % normalize

bases = [zeta^2*w_ratio^2, ...
         zeta^2*w_ratio^2 , ...
         (1 + w_ratio^4 + 4*zeta^2*w_ratio^2 - 2*w_ratio^2), ...
         (w_ratio^2 - 1), ...
         1];

for abc=1:length(quad)
    matrix_tmp = coeffs(quad(abc),[r_k r_b],'all');
    matrix(abc,:) = matrix_tmp([7 4 1 2 3]);
end
matrix = matrix ./ bases;
matrix = simplify(matrix,'IgnoreAnalyticConstraints',true)

latex(matrix)
bases_rk = bases .* [r_k^2, r_k^2*r_b, r_k^2*r_b^2, r_k*r_b^2, r_b^2];
pretty(bases_rk')
latex(bases_rk)

% verify that setting m=1 (so a*1^2 + b*1 + c = a+b+c = 0) causes f=1
eqn_m_1 = sum(quad) == 0;
soln_f = solve(eqn_m_1,f);
assert(soln_f == 1)

%% optimal powertrain case
bases_rk = subs(bases_rk,{w_ratio,r_b},{1,2}).';

abc = matrix * bases_rk/4;
pretty(abc)

discrim = sqrt(abc(2)^2 - 4 * abc(1) * abc(3));
m_plus  = simplify((-abc(2) + discrim) / (2*abc(1)));
m_minus = simplify((-abc(2) - discrim) / (2*abc(1)));

syms r_k_zeta
m_plus  = subs(m_plus, r_k, r_k_zeta/zeta);
m_minus = subs(m_minus,r_k, r_k_zeta/zeta);

pretty(m_plus)
pretty(m_minus)

%% sub in numerical
% relationship between m, f_sat, and w for the ideal case (max power control)
% see notebook p115-118
close all

f_sat = linspace(0.01,1);
r_k_z = linspace(0.01,3);

[F,RKZ] = meshgrid(f_sat,r_k_z);

mult_num = vpa(subs(m_plus,{f,r_k_zeta},{F,RKZ}));
mult_num(:,:,2) = vpa(subs(m_minus,{f,r_k_zeta},{F,RKZ}));

energy = F.^2 ./ mult_num(:,:,1);

figure
contourf(F,RKZ,mult_num(:,:,1))
xlabel('f_{sat}')
ylabel('|r_k \zeta|')
title('m_+')
colorbar

figure
contourf(F,RKZ,mult_num(:,:,2),[min(mult_num(:,:,2),[],'all') 0:.1:1 max(mult_num(:,:,2),[],'all')])
xlabel('f_{sat}')
ylabel('|r_k \zeta|')
title('m_-')
colorbar
caxis([0 1])

figure
contourf(F,RKZ,energy)
xlabel('f_{sat}')
ylabel('|r_k \zeta|')
title('Energy scale for m_-')
colorbar
caxis([0 1])

r_k_z_selected = r_k_z([1 4 7 11 17 34 100]); % corresponds to values [.01 .1 .2 .3 .5 1 3]
figure;
hold on

colors = {'m','r','y','g','c','b','k'};

for i=1:length(r_k_z_selected)
    rkzi = r_k_z_selected(i);
    idx = RKZ == rkzi;
    f_i = F(idx);
    m_tmp = mult_num(:,:,1);
    m_i = m_tmp(idx);
    e_i = energy(idx);

    subplot 121
    plot(f_i,m_i,colors{i},'DisplayName',['|r_k \zeta| = ' num2str(rkzi)])
    hold on
    xlabel('force saturation factor f_{sat}')
    ylabel('powertrain control multiplier m')

    subplot 122
    plot(f_i,e_i,colors{i},'DisplayName',['|r_k \zeta| = ' num2str(round(rkzi,1))])
    hold on
    xlabel('force saturation factor f_{sat}')
    ylabel('energy multiplier e')
end
legend
sgtitle('Effect of Hydrodynamics on Force Saturation Behavior')
improvePlot
