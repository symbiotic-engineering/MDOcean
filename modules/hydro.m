function [F_heave, F_surge, m, D_env] = hydro(t, s, m_float, D_sft, rho_w, g, Hs, T, t_f)
% unused: d_WEC, N_WEC, d_farm, d_shore

w = 2*pi/T;         % frequency
k = w^2 / g;        % wave number
V_g = g /(2*w);     % group velocity

r = D_sft / 2;      % radius
A_w = pi * r^2;     % waterplane area
draft = t_f / 2;    % estimated depth
vol = A_w * draft;  % submerged volume

% hydrodynamic coefficients
vec_size = ones(1,size(s,2));
A       = 1/2 * rho_w * 4/3 * pi * r^3 * 0.63 * vec_size; % added mass
gamma   = rho_w * g * A_w * vec_size; % Froude Krylov / diffraction
B       = k / (4 * rho_w * g * V_g) * gamma.^2; % radiation damping

wave = Hs * sin(w * t); % in reality this is a function of d_shore

F_hydrostatic = -rho_w * g * A_w * s(1,:);
F_radiation = -B .* s(2,:);
F_excitation = gamma .* wave;
F_heave = F_hydrostatic + F_radiation + F_excitation;

F_surge = 2 * k * rho_w * g * vol * (1 - exp(k*draft));

m = m_float + A;

D_env = 0.5 * vec_size; % in reality this is a function of d_farm and forces

end

