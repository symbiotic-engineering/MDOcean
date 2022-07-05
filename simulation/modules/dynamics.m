
function [F_heave_max, F_surge_max, F_ptrain_max, P_var, P_elec, P_matrix, h_s_extra] = dynamics(in,m_float,V_d,draft)

% use probabilistic sea states for power
[Hs,T] = meshgrid(in.T,in.Hs);
P_matrix = get_power_force(in,T,Hs,m_float,V_d,draft);

% weight power across all sea states
P_weighted = P_matrix .* in.JPD / 100;
P_elec = sum(P_weighted(:)); 

% use max sea states for structural forces and max amplitude
[~,F_heave_max,F_surge_max,...
    F_ptrain_max,h_s_extra] = get_power_force(in, ...
                            in.T_struct, in.Hs_struct, m_float, V_d, draft);

% coefficient of variance (normalized standard deviation) of power
P_var = std(P_matrix(:), in.JPD(:)) / P_elec;
P_var = P_var * 100; % convert to percentage

end

function [P_matrix, F_heave, F_surge, F_ptrain, F_ptrain_max, h_s_extra] = get_power_force(in,T,Hs, m_float,V_d, draft)
    % get unsaturated response
    [w,A,B_h,K_h,Fd,k_wvn] = dynamics_simple(Hs, T, in.D_f, in.rho_w, in.g);
    m = m_float + A;
    b = B_h + in.B_p;
    k = in.w_n^2 * m;
    K_p = k - K_h;
    X_unsat = get_response(w,m,b,k,Fd);
    F_ptrain_over_x = sqrt( (in.B_p * w).^2 + (K_p).^2 );
    F_ptrain_unsat = F_ptrain_over_x .* X_unsat;
    
    % get saturated response
    r = min(in.F_max ./ F_ptrain_unsat, 1);%fcn2optimexpr(@min, in.F_max ./ F_ptrain_unsat, 1);
    alpha = 2/pi * ( 1./r .* asin(r) + sqrt(1 - r.^2) );
    f_sat = alpha .* r;
    mult = get_multiplier(f_sat,m,b,k,w, B_h/in.B_p, K_h/K_p);
    b_sat = B_h + mult * in.B_p;
    k_sat = K_h + mult * K_p;
    X_sat = get_response(w,m,b_sat,k_sat,Fd);
    
    % calculate power
    P_matrix = 1/2 * (mult * in.B_p) .* w.^2 .* X_sat.^2;
    
    % calculate forces
    if nargout > 1
        F_ptrain = mult .* F_ptrain_over_x .* X_sat;
        F_ptrain_max = max(F_ptrain,[],'all');
        F_err_1 = abs(F_ptrain ./ (in.F_max * alpha) - 1);
        F_err_2 = abs(F_ptrain ./ (f_sat * F_ptrain_unsat) - 1);
        % 20 percent error acceptable for now
        if any(f_sat<1,'all')
            assert(F_err_1(f_sat < 1) < 0.20);
        end
        assert(F_err_2 < 0.20);

        F_heave_fund = sqrt( (mult * in.B_p * w).^2 + (mult * K_p - m_float * w.^2).^2 ) .* X_sat; % includes powertrain force and D'Alembert force
        F_heave = min(F_heave_fund, in.F_max + m_float * w.^2 .* X_sat);
        %assert(F_heave <= in.F_max);

        F_surge = Hs * in.rho_w * in.g * V_d .* (1 - exp(-k_wvn*draft));
        X_max = max(X_sat,[],'all');
        h_s_extra = (in.h_s - in.T_s - (in.h_f - in.T_f) - X_max) / in.h_s; % extra height on spar after accommodating float displacement
    end
end

function [w,A,B,K,Fd,k] = dynamics_simple(Hs, T, D_f, rho_w, g)
    w = 2*pi./T;         % frequency
    k = w.^2 / g;        % wave number
    V_g = g ./(2*w);     % group velocity

    r = D_f / 2;      % radius
    A_w = pi * r^2;     % waterplane area
    
    A       = 1/2 * rho_w * 4/3 * pi * r^3 * 0.63; % added mass
    gamma   = rho_w * g * A_w; % Froude Krylov / diffraction
    B       = k ./ (4 * rho_w * g * V_g) * gamma.^2; % radiation damping
    K       = rho_w * g * A_w;  % hydrostatic stiffness
    Fd      = gamma * Hs;       % excitation force of wave
end

function X = get_response(w,m,b,k,Fd)
    real_term = b.*w;
    imag_term = k - m*w.^2;
    X_over_F_mag = ((real_term).^2 + (imag_term).^2).^(-1/2);
    %X_over_F_phase = atan2(imag_term,real_term);
    X = X_over_F_mag .* Fd;
end

function mult = get_multiplier(f_sat,m,b,k,w,r_b,r_k)
    % m, k, and r_k are scalars.
    % All other inputs are 2D arrays, the dimension of the sea state matrix.
    % This algebra is written out on p198 of my notebook.

    % speedup: only do math for saturated sea states, since unsat will = 1
    idx_no_sat = f_sat == 1;
    f_sat(idx_no_sat) = NaN;
    b(idx_no_sat) = NaN;
    w(idx_no_sat) = NaN;
    r_b(idx_no_sat) = NaN;
    
    % coefficients defined for convenience
    alpha_b = 1./(r_b + 1);
    alpha_k = 1./(r_k + 1);
    beta_b = r_b./(r_b + 1);
    beta_k = r_k./(r_k + 1);
    
    m2_w4 = m.^2 .* w.^4;
    two_k_m_w2 = 2 * k .* m .* w.^2;
    inv_f_sat2 = f_sat.^(-2);
    bw2 = (b.*w).^2;
    k2 = k.^2;

    % quadratic formula coeffs: 0 = a_quad * mult^2 + b_quad * mult + c_quad
    a_quad = (alpha_b.^2 - inv_f_sat2) .* bw2 ...
           + (alpha_k.^2 - inv_f_sat2) .* k2 ...
           +  two_k_m_w2 - m2_w4;
    b_quad = 2 * (alpha_b .* beta_b .* bw2 ...
                +  alpha_k .* beta_k .* k2) ...
                -  alpha_k .* two_k_m_w2;
    c_quad = beta_b.^2 .* bw2 + beta_k.^2 .* k2 ...
           - beta_k .* two_k_m_w2 + m2_w4;

    % solve the quadratic formula
    determinant = sqrt(b_quad .^ 2 - 4 * a_quad .* c_quad);
    num = -b_quad + determinant;
    num(:,:,2) = -b_quad - determinant;
    den = 2 * a_quad;
    roots = num ./ den;

    % choose which of the two roots to use
    which_soln = roots == real(roots) & roots > 0 & roots <= 1; % real solns on (0, 1]
    both_ok = sum(which_soln,3) == 2;
    
    which_soln(idx_no_sat) = 1; % temporarily mark the non-saturated solutions
                                % as having one solution, to ensure the 
                                % logic below works correctly

    if any(both_ok,'all')
        warning('Using outliers to determine relevant quadratic formula solution')
        [row,col]=find(both_ok);

        % if both are ok, choose the first one arbitrarily for now
        which_soln(row,col,2) = false; 
        mult_1 = get_relevant_soln(which_soln,roots,idx_no_sat);   

        % now choose the second one arbitrarily
        which_soln(row,col,2) = true;
        which_soln(row,col,1) = false;
        mult_2 = get_relevant_soln(which_soln,roots,idx_no_sat);   

        % compare outliers to figure out if first or second is better
        window = 6;
        outliers_1 = isoutlier(mult_1,'movmedian',window);
        outliers_2 = isoutlier(mult_2,'movmedian',window);
        
        outliers_1_relevant = outliers_1(both_ok);
        outliers_2_relevant = outliers_2(both_ok);

        one_all_outliers = all(outliers_1_relevant);
        two_all_outliers = all(outliers_2_relevant);
        one_all_ok = all(~outliers_1_relevant);
        two_all_ok = all(~outliers_2_relevant);
        
        use_1 = (two_all_outliers && ~one_all_outliers) || (one_all_ok && ~two_all_ok);
        use_2 = (one_all_outliers && ~two_all_outliers) || (two_all_ok && ~one_all_ok);
%         use_1 = sum(outliers_1_relevant,'all') < sum(outliers_2_relevant,'all');
%         use_2 = sum(outliers_1_relevant,'all') > sum(outliers_2_relevant,'all');

        if use_1
            mult = mult_1;
        elseif use_2
            mult = mult_2;
        else
            figure
            subplot 121
            contourf(mult_1)
            subplot 122
            contourf(mult_2)

            error(['Failed to figure out which solution to the quadratic' ...
                'equation is relevant, try manual inspection.'])
        end
        
    else
        
        num_solns = sum(which_soln,3);
        assert(all( num_solns == 1,'all') ); % confirm that 1 soln per sea state meets criteria

        mult = get_relevant_soln(which_soln,roots,idx_no_sat);   

    end

    assert(all(~isnan(mult),'all'))
end

function mult = get_relevant_soln(which_soln, roots, idx_no_sat)
% pick the specified roots using multidimensional logical indexing

    mult = NaN(size(idx_no_sat));

    % figure out 3d and 2d indices
    idx_3d_first_sol = which_soln;
    idx_3d_first_sol(:,:,2) = false;
    idx_3d_second_sol = which_soln;
    idx_3d_second_sol(:,:,1) = false;
    idx_2d_first_sol = which_soln(:,:,1);
    idx_2d_second_sol = which_soln(:,:,2);
    
    mult(idx_2d_first_sol) = roots(idx_3d_first_sol);
    mult(idx_2d_second_sol) = roots(idx_3d_second_sol);
    mult(idx_no_sat) = 1;
end
