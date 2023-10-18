clear
close all

% sweep variables
zeta_u = 1;
w_n_u = 1;
X_prime = 1;
w = 1;
F_max_over_Fp = 1;

% constant variables - changing this should have no effect
m_set = 1;

% parameters
p = struct('Kh',w_n_u^2, 'Bh',2*zeta_u*w_n_u, 'm',m_set, 'w',w, ...
           'Fh',X_prime*m_set*w_n_u^2, 'Fmax',.4);

% impedance matching
w_star = 1;
r_b = 2;
p.Kp = p.m * p.w^2 - p.Kh; % this should techncially depend on w*
p.Bp = p.Bh; % this should techncially depend on rb

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

% equation of motion
function [ydot,Fp,P] = dynamics(t,y,p)

% states
x = y(1,:);
xdot = y(2,:);

% forces
Fh = p.Fh * sin(p.w * t);
Fhyd = p.Kh * x + p.Bh * xdot;
Fp = p.Kp * x + p.Bp * xdot;
Fp = min(max(Fp/p.Fmax, -1), 1) * p.Fmax;
F = Fh - Fp - Fhyd;

% state derivatives
xddot = F / p.m;
ydot = [xdot; xddot];

% power
P = Fp .* xdot;

end