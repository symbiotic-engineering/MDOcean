p = parameters();

p.Hs  = 1 + 0 * p.Hs;
p.JPD = 1 + 0 * p.JPD;
p.T   = 8 + 0 * p.T;
p.use_force_sat = false;

b = var_bounds();
X = [b.X_noms; 1];
[~, ~, ~, val] = simulation(X,p);

K_mult = linspace(.8,1.2,15);
B_mult = linspace(0.5,1.5,14).';

% copy paste these lines at the end of controller in dynamics.m
% B_p_0 = B_p(1);
% K_p_0 = K_p(1);
% B_mult = linspace(.5,1.5,size(w,1)).';
% K_mult = linspace(.8,1.2,size(w,2));
% B_p = repmat(B_mult,[1 size(w,2)]) * B_p_0; 
% K_p = repmat(K_mult,[size(w,1) 1]) * K_p_0;

figure
contourf(K_mult,B_mult,val.P_mech/p.power_scale_multibody)
xlabel('K_p multiplier');
ylabel('B_p multiplier');
title('Power for various control');
colorbar
grid on

hold on
plot(1,1,'ro')