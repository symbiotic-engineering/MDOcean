function spar_exc = get_spar_exc(g)

hydro = struct();
hydro = readWAMIT(hydro,'rm3.out',[]); % function from WECSim

%B_over_rho_w = hydro.B(9,9,:);
%B_over_rho_w = B_over_rho_w(:);

gamma_over_rho_g = hydro.ex_ma(9,1,:);
gamma_over_rho_g = gamma_over_rho_g(:);

gamma_phase= -hydro.ex_ph(9,1,:);
gamma_phase = gamma_phase(:);

k = hydro.w.^2 / g;

spar_exc = struct('gamma_over_rho_g',gamma_over_rho_g,...
                  'gamma_phase',gamma_phase,...
                  'k', k, ...
                  'T_s',   35 ); % spar draft corresponding to the spar excitation coeffs
            