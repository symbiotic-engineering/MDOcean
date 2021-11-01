
function [F_heave, F_surge, F_ptrain, D_env, P_elec, P_matrix, F_ptrain_unsat] = dynamicSimulation(x,p,m_float,t_f)
   
[Hs,T] = meshgrid(p.T,p.Hs);

% get unsaturated response
[w,A,B,K,Fd] = dynamics_simple(Hs, T, x.D_sft,  p.rho_w, p.g);       
m = m_float + A;
b = B + x.D_int;
k = x.w_n^2 * m;
K_int = k - K;
[~,F_ptrain_unsat] = get_response(w,m,b,k,Fd);

mult = min(p.F_max ./ F_ptrain_unsat, 1);
% get saturated response
% fixme: should multiply mult (saturation multiplier) by a fourier multiplier to get total mult
b_sat = B + mult * x.D_int;
k_sat = K + mult * K_int;
[X_sat,F_ptrain] = get_response(w,m,b_sat,k_sat,Fd);

P_matrix = 1/2 * (mult * x.D_int) .* w.^2 .* X_sat.^2;

% weight power across all sea states
P_weighted = P_matrix .* p.JPD;
P_elec = mean(P_weighted(:));

% covert time series to scalar outputs
D_env = 0; 
F_heave = 10;%max(F_heave);
F_surge = 10;% fixme this should use structural sea state
F_ptrain = mean(F_ptrain(:)); % fixme this should use structural sea state

end

function [w,A,B,K,Fd] = dynamics_simple(Hs,T,D_sft, rho_w, g)
    w = 2*pi./T;         % frequency
    k = w.^2 / g;        % wave number
    V_g = g ./(2*w);     % group velocity

    r = D_sft / 2;      % radius
    A_w = pi * r^2;     % waterplane area
    
    A       = 1/2 * rho_w * 4/3 * pi * r^3 * 0.63; % added mass
    gamma   = rho_w * g * A_w; % Froude Krylov / diffraction
    B       = k ./ (4 * rho_w * g * V_g) * gamma.^2; % radiation damping
    K       = rho_w * g * A_w;  % hydrostatic stiffness
    Fd      = gamma * Hs;       % excitation force of wave
end

function [X,F_ptrain] = get_response(w,m,b,k,Fd)
    real_term = b.*w/m;
    imag_term = k/m - w.^2;
    X_over_F_mag = ((real_term).^2 + (imag_term).^2).^(-1/2);
    %X_over_F_phase = atan2(imag_term,real_term);
    X = X_over_F_mag .* Fd;
    F_ptrain = (b.*w+k).*X;
end
