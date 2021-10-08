function [F_hydro, m, D_env] = hydro(t, s, m_float, D_sft, d_WEC, N_WEC, d_farm, d_shore, rho_w, g, Hs, T)

w = 2*pi/T;     % frequency
r = D_sft / 2;  % radius
A_w = pi * r^2; % waterplane area
k = w^2 / g;    % wave number
V_g = g /(2*w); % group velocity

% hydrodynamic coefficients
A       = 1/2 * rho_w * 4/3 * pi * r^3 * 0.63 * ones(1,size(s,2)); % added mass
gamma   = rho_w * g * A_w * ones(1,size(s,2)); % Froude Krylov / diffraction
B       = k / (4 * rho_w * g * V_g) * gamma.^2; % radiation damping

wave = Hs * sin(w * t); % in reality this is a function of d_shore

F_hydrostatic = -rho_w * g * A_w * s(1,:);
F_radiation = -B .* s(2,:);
F_excitation = gamma .* wave;
F_hydro = F_hydrostatic + F_radiation + F_excitation;

m = m_float + A;

D_env = 0.5 * ones(1,size(s,2)); % in reality this is a function of d_farm and forces

end

