clear
close all

% sweep variables
zeta_u = [.2 .5 .8 1]; % sets level of damping of uncontrolled system
w_u_star = 0.5 : 0.1 : 0.9; % sets frequency ratio of uncontrolled system
F_max_over_Fp = 0.5 : 0.1 : 1.2; % sets level of saturation

% generate sweep grid
[ZETA_U,W_U_STAR,F_MAX_FP] = meshgrid(zeta_u, w_u_star, F_max_over_Fp);

% constant variables - changing these should have no effect on the nondimensional outputs
m = 1;  % dimensional term for mass
w = 1;  % dimensional term for frequency
F_h = 1; % dimensional term for exciting force

% get results
[avg_pwr, max_x, max_xdot, pwr_ratio, x_ratio, xdot_ratio] = describing_fcn(ZETA_U, W_U_STAR, F_MAX_FP, m, w, F_h);
[avg_pwr2, max_x2, max_xdot2, pwr_ratio2, x_ratio2, xdot_ratio2] = ground_truth(  ZETA_U, W_U_STAR, F_MAX_FP, m, w, F_h);

% combine results
dim_cat = ndims(ZETA_U)+1;
results_sim = cat(dim_cat, avg_pwr,  max_x,  max_xdot,  pwr_ratio,  x_ratio,  xdot_ratio);
results_act = cat(dim_cat, avg_pwr2, max_x2, max_xdot2, pwr_ratio2, x_ratio2, xdot_ratio2);
results_str = {'Average Power','Max X','Max Xdot','Power Ratio','X Ratio','Xdot Ratio'};

results_unsat_sim = results_sim(F_MAX_FP >= 1);
results_unsat_act = results_act(F_MAX_FP >= 1);

for j = 1:size(results_sim,dim_cat)
    figure
    res_sim = results_sim(:,:,:,j);
    res_act = results_act(:,:,:,j);
    color_lims = [min([res_sim,res_act],[],'all') max([res_sim,res_act],[],'all')];

    % plot simulated results
    subplot 131
    slice(ZETA_U,W_U_STAR,F_MAX_FP, res_sim,...
        zeta_u(2:end),w_u_star(end),F_max_over_Fp,'nearest');
    caxis(color_lims)
    colorbar
    xlabel('\zeta_u')
    ylabel('\omega_u^*')
    zlabel('F_{max}/F_p')
    title('Describing Function')

    % plot actual results
    subplot 132
    slice(ZETA_U,W_U_STAR,F_MAX_FP, res_act,...
        zeta_u(2:end),w_u_star(end),F_max_over_Fp,'nearest');
    caxis(color_lims)
    colorbar
    xlabel('\zeta_u')
    ylabel('\omega_u^*')
    zlabel('F_{max}/F_p')
    title('Numerical Integration')

    error_h = subplot(1,3,3);
    slice(ZETA_U,W_U_STAR,F_MAX_FP, (res_sim-res_act)./res_act,...
        zeta_u(2:end),w_u_star(end),F_max_over_Fp,'nearest');
    colormap(error_h,bluewhitered)
    colorbar
    xlabel('\zeta_u')
    ylabel('\omega_u^*')
    zlabel('F_{max}/F_p')
    title('Fractional Error')

    sgtitle(results_str{j})
end

function [avg_pwr, max_x, max_xdot, ...
    pwr_ratio, x_ratio, xdot_ratio] = describing_fcn(zeta_u, w_u_star, ...
                                                    F_max_over_Fp, m, w, F_h)

r_b = 2;
w_star = 1;

F_max_over_Fp = min(F_max_over_Fp, 1);
f_sat = 2/pi * (F_max_over_Fp .* sqrt(1 - F_max_over_Fp.^2) + asin(F_max_over_Fp));
rkz = abs(w_u_star / w_star * r_b .* zeta_u ./ ( ((w_u_star / w_star).^2 - 1) * (r_b-1)));
% m_sat equation assumes r_b = 2 and w_star = 1
m_sat = -f_sat .* (f_sat .* rkz.^2 - f_sat + 2*rkz .* sqrt(-f_sat.^2 + rkz.^2 + 1)) ./ (f_sat.^2 .* rkz.^2 + f_sat.^2 - 4*rkz.^2);
e = f_sat.^2 ./ m_sat;

% P_unsat equation assumes r_b = 2 and w_star = 1
P_unsat = 1/16 * F_h^2 / (m*w) * w_u_star ./ zeta_u;
P_sat = P_unsat .* e;

zeta = r_b/(r_b - 1) * w_star ./ w_u_star .* zeta_u;
X_unsat = F_h^2 / (m*w^2) * w_star.^2 ./ sqrt( (1 - w_star^2)^2 + (2*zeta*w_star).^2);
X_sat = X_unsat .* f_sat ./ m_sat;

avg_pwr = P_sat;
max_x = X_sat;
max_xdot = max_x * w;
pwr_ratio = e;
x_ratio = f_sat ./ m_sat;
xdot_ratio = f_sat ./ m_sat;

end

function [avg_pwr, max_x, max_xdot, pwr_ratio, x_ratio, xdot_ratio] = ground_truth(ZETA_U, W_U_STAR, F_MAX_FP, m, w, F_h)

[avg_pwr, max_x, max_xdot, pwr_ratio, x_ratio, xdot_ratio] = deal(zeros(size(ZETA_U)));

for i = 1:numel(ZETA_U)
    zeta_u = ZETA_U(i);
    w_u_star = W_U_STAR(i);
    F_max_over_Fp = F_MAX_FP(i);

    [avg_pwr(i), max_x(i), max_xdot(i)] = run_ode(zeta_u, w_u_star, F_max_over_Fp, m, w, F_h, false);
    [avg_pwr_unsat, max_x_unsat, max_xdot_unsat] = run_ode(zeta_u, w_u_star, 1e8, m, w, F_h, false);
    
    pwr_ratio(i) = avg_pwr(i) / avg_pwr_unsat;
    x_ratio(i) = max_x(i) / max_x_unsat;
    xdot_ratio(i) = max_xdot(i) / max_xdot_unsat;
end

end

function [avg_pwr, max_x, max_xdot] = run_ode(zeta_u, w_u_star, F_max_over_Fp, m, w, F_h, plotOn)

% dependent variables
w_n_u = w ./ w_u_star;

% parameters
p = struct('Kh',w_n_u^2*m, 'Bh',2*zeta_u*w_n_u*m, 'm',m, 'w',w, ...
           'Fh',F_h);

% impedance matching
w_star = 1;
r_b = 2;
p.Kp = p.m * p.w^2 - p.Kh; % fixme: this should techncially depend on w*
p.Bp = p.Bh; % fixme: this should techncially depend on rb

% maximum force
wn = p.w / w_star;
w_u_star = p.w / w_n_u;
zeta = r_b / (r_b - 1) * w_star / w_u_star * zeta_u;
X = p.Fh / (p.m * wn^2) / sqrt( (1 - w_star^2)^2 + (2*zeta*w_star)^2 );
F_p = ((p.Bp * p.w)^2 + p.Kp^2) * X;
p.F_max = F_max_over_Fp * F_p;

% ode inputs
T = 2*pi/p.w;
y0 = [0,0];
tspan = linspace(0,5*T,501);

% ode solve
[t,y] = ode45(@(t,y)dynamics(t,y,p),tspan,y0);
[~,Fp,P] = dynamics(t',y',p);

% plot
if plotOn
    figure
    plot(t,y)
    hold on
    plot(t,Fp,t,P)
    legend('x','xdot','Fp','P')
end

avg_pwr = mean(P(t>=4*T));
max_x = max(abs(y(t>=4*T,1)));
max_xdot = max(abs(y(t>=4*T,2)));

end 

% equation of motion
function [ydot,Fp,P] = dynamics(t,y,p)

% states
x = y(1,:);
xdot = y(2,:);

% forces
Fh = p.Fh * sin(p.w * t);
Fhyd = p.Kh * x + p.Bh * xdot;
Fp = p.Kp * x + p.Bp * xdot;
Fp = min(max(Fp/p.F_max, -1), 1) * p.F_max;
F = Fh - Fp - Fhyd;

% state derivatives
xddot = F / p.m;
ydot = [xdot; xddot];

% power
P = Fp .* xdot;

end