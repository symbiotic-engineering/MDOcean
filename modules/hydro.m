function [F_hydro, m, D_env] = hydro(t, s, L, W, M, d_WEC, N_WEC, d_farm, d_shore, rho_w, g, Hs, T)

% hydrodynamic coefficients (in reality, these are functions of s, L, W, d_WEC, N_WEC)
B = 1; % radiation damping
A = 1; % added mass
gamma = 1; % diffraction

wave = Hs * sin(T * t); % in reality this is a function of d_shore

F_hydrostatic = -rho_w * g * L * W * s(1);
F_radiation = -B * s(2);
F_excitation = gamma * wave;
F_hydro = F_hydrostatic + F_radiation + F_excitation;

rho = 1000; % in reality this is a function of material M
mass = L * W * rho;
m = mass + A;

D_env = 0.5; % in reality this is a function of d_farm and forces

end

