clear
close all

% sweep variables
zeta_u = 1; % sets level of damping of uncontrolled system
w_u_star = .5; % sets frequency ratio of uncontrolled system
F_max_over_Fp = 1; % sets level of saturation

% constant variables - changing these should have no effect on the nondimensional outputs
m = 1;  % dimensional term for mass
w = 1;  % dimensional term for frequency
F_h = 1; % dimensional term for exciting force

[avg_pwr, max_x, max_xdot, pwr_ratio, x_ratio, xdot_ratio] = describing_fcn(zeta_u, w_u_star, F_max_over_Fp, m, w, F_h)
[avg_pwr, max_x, max_xdot, pwr_ratio, x_ratio, xdot_ratio] = ground_truth(  zeta_u, w_u_star, F_max_over_Fp, m, w, F_h)

function [avg_pwr, max_x, max_xdot, pwr_ratio, x_ratio, xdot_ratio] = describing_fcn(zeta_u, w_u_star, F_max_over_Fp, m, w, F_h)

r_b = 2;
w_star = 1;

f_sat = 2/pi * (F_max_over_Fp * sqrt(1 - F_max_over_Fp^2) + asin(F_max_over_Fp));
rkz = abs(w_u_star / w_star * r_b * zeta_u / ( ((w_u_star / w_star)^2 - 1) * (r_b-1)));
% m_sat equation assumes r_b = 2 and w_star = 1
m_sat = -f_sat * (f_sat * rkz^2 - f_sat + 2*rkz * sqrt(-f_sat^2 + rkz^2 + 1)) / (f_sat^2 * rkz^2 + f_sat^2 - 4*rkz^2);
e = f_sat^2 / m_sat;

% P_unsat equation assumes r_b = 2 and w_star = 1
P_unsat = 1/16 * F_h^2 / (m*w) * w_u_star / zeta_u;
P_sat = P_unsat * e;

zeta = r_b/(r_b - 1) * w_star / w_u_star * zeta_u;
X_unsat = F_h^2 / (m*w^2) * w_star^2 / sqrt( (1 - w_star^2)^2 + (2*zeta*w_star)^2);
X_sat = X_unsat * f_sat / m_sat;

avg_pwr = P_sat;
max_x = X_sat;
max_xdot = max_x * w;
pwr_ratio = e;
x_ratio = f_sat / m_sat;
xdot_ratio = f_sat / m_sat;

end

function [avg_pwr, max_x, max_xdot, pwr_ratio, x_ratio, xdot_ratio] = ground_truth(zeta_u, w_u_star, F_max_over_Fp, m, w, F_h)

[avg_pwr, max_x, max_xdot] = run_ode(zeta_u, w_u_star, F_max_over_Fp, m, w, F_h);
[avg_pwr_unsat, max_x_unsat, max_xdot_unsat] = run_ode(zeta_u, w_u_star, 1e8, m, w, F_h);

pwr_ratio = avg_pwr / avg_pwr_unsat;
x_ratio = max_x / max_x_unsat;
xdot_ratio = max_xdot / max_xdot_unsat;

end

function [avg_pwr, max_x, max_xdot] = run_ode(zeta_u, w_u_star, F_max_over_Fp, m, w, F_h)

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
tspan = [0 5*T];

% ode solve
[t,y] = ode45(@(t,y)dynamics(t,y,p),tspan,y0);
[~,Fp,P] = dynamics(t',y',p);

% plot
figure
plot(t,y)
hold on
plot(t,Fp,t,P)
legend('x','xdot','Fp','P')

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