% symbolic gradient for power

syms P C w T k V_g g r D A_w Fd m K rho_w Hs A gamma B TF real
syms s

w = 2*pi/T;         % frequency
k = w^2 / g;        % wave number
V_g = g /(2*w);     % group velocity

r = D / 2;          % radius
A_w = pi * r^2;     % waterplane area

A       = 1/2 * rho_w * 4/3 * pi * r^3 * 0.63; % added mass
gamma   = rho_w * g * A_w; % Froude Krylov / diffraction
B       = k / (4 * rho_w * g * V_g) * gamma.^2; % radiation damping

s = 1i * w;
TF = 1/abs( s^2*(m+A) + s*(B+C) + rho_w*g*A_w + K );

Fd = Hs * gamma;

P = simplify(1/2 * C * w^2 * Fd^2 * TF^2)

gradP = gradient(P, [C T Hs D m K])

dir = ['generated',filesep,'objective.m'];

matlabFunction(P,gradP,'file',dir,'vars',[C T Hs D m K rho_w g]);
