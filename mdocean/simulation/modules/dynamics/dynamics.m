
function [F_heave_max, F_surge_max, F_ptrain_max, P_var, P_elec, P_matrix, h_s_extra] = dynamics(in,m_float,V_d,draft)

    % use probabilistic sea states for power
    [T,Hs] = meshgrid(in.T,in.Hs);
    [P_matrix,h_s_extra] = get_power_force(in,T,Hs,m_float,V_d,draft);
    
    % account for powertrain electrical losses
    P_matrix = P_matrix * in.eff_pto;
    
    % saturate maximum power
    P_matrix = min(P_matrix,in.power_max);
    
    % weight power across all sea states
    P_weighted = P_matrix .* in.JPD / 100;
    P_elec = sum(P_weighted(:)); 
    
    assert(isreal(P_elec))
    
    % use max sea states for structural forces and max amplitude
    [~,~,F_heave_max,F_surge_max,...
        F_ptrain_max] = get_power_force(in, ...
                                in.T_struct, in.Hs_struct, m_float, V_d, draft);
    
    % coefficient of variance (normalized standard deviation) of power
    P_var = std(P_matrix(:), in.JPD(:)) / P_elec;
    P_var = P_var * 100; % convert to percentage

end

function [P_matrix, h_s_extra, F_heave, F_surge, F_ptrain_max] = get_power_force(in,T,Hs, m_float,V_d, draft)
    % get unsaturated response
    [w,A,B_h,K_h,Fd,k_wvn] = dynamics_simple(Hs, T, in.D_f, in.T_f, in.rho_w, in.g);
    m = m_float + A;
    b = B_h + in.B_p;
    k = in.w_n^2 * m;
    K_p = k - K_h;
    X_unsat = get_response(w,m,b,k,Fd);

    % confirm unsaturated response doesn't exceed maximum capture width
%     P_unsat = 1/2 * in.B_p * w.^2 .* X_unsat.^2;
%     P_wave = in.rho_w * in.g^2 / (64*pi) * T .* Hs.^2;
%     CW = P_unsat ./ P_wave;
%     CW_max = in.g * T.^2 / (4*pi^2);
%     assert( all(CW <= CW_max, 'all') ); % check hydro coeffs if this fails

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
    
    X_max = max(X_sat,[],'all');
    h_s_extra = (in.h_s - in.T_s - (in.h_f - in.T_f) - X_max) / in.h_s; % extra height on spar after accommodating float displacement

    % calculate forces
    if nargout > 2
        F_ptrain = mult .* F_ptrain_over_x .* X_sat;
        F_ptrain_max = max(F_ptrain,[],'all');
        F_err_1 = abs(F_ptrain ./ (in.F_max * alpha) - 1);
        F_err_2 = abs(F_ptrain ./ (f_sat * F_ptrain_unsat) - 1);
        % 0.1 percent error
        if any(f_sat<1,'all')
            assert(F_err_1(f_sat < 1) < 1e-3);
        end
        assert(F_err_2 < 1e-3);

        F_heave_fund = sqrt( (mult * in.B_p * w).^2 + (mult * K_p - m_float * w.^2).^2 ) .* X_sat; % includes powertrain force and D'Alembert force
        F_heave = min(F_heave_fund, in.F_max + m_float * w.^2 .* X_sat);
        %assert(F_heave <= in.F_max);

        F_surge = Hs * in.rho_w * in.g * V_d .* (1 - exp(-k_wvn*draft));
    end
end

function [w,A,B,K,Fd,k] = dynamics_simple(Hs, T, D_f, T_f, rho_w, g)
    w = 2*pi./T;        % angular frequency
    k = w.^2 / g;       % wave number (dispersion relation for deep water)

    r = D_f / 2;        % radius
    draft = T_f;        % draft below waterline
    A_w = pi * r^2;     % waterplane area
    
    [A_over_rho, B_over_rho_w, gamma_over_rho_g] = get_hydro_coeffs(r, k, draft);  
    
    A = rho_w * A_over_rho;                 % added mass
    B = rho_w * w .* B_over_rho_w;          % radiation damping
    gamma   = rho_w * g * gamma_over_rho_g; % froude krylov coefficient
    K       = rho_w * g * A_w;              % hydrostatic stiffness
    Fd      = gamma .* Hs;                  % excitation force of wave
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

    % speedup: only do math for saturated sea states, since unsat will = 1
    idx_no_sat = f_sat == 1;
    f_sat(idx_no_sat) = NaN;
    b(idx_no_sat) = NaN;
    w(idx_no_sat) = NaN;
    r_b(idx_no_sat) = NaN;
    
    [a_quad, b_quad, c_quad]  = get_abc_symbolic(f_sat,m,b,k,w,r_b,r_k);

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
        
        mult = handle_two_solns(both_ok,which_soln,roots,idx_no_sat,a_quad,b_quad,c_quad);
    else    
        num_solns = sum(which_soln,3);
        if ~(all( num_solns == 1,'all') )
            which_soln(num_solns==0) = roots(num_solns==0) > 0 & roots(num_solns==0) <= 1.001; % wider tolerance
            num_solns(num_solns==0) = sum(which_soln(num_solns==0),3);
            if ~(all( num_solns == 1,'all'))
                % if still a problem, proceed with which_soln set to zero
                % for the problem sea states, which will set mult = 0
                disp('ohno')
            end
            % confirm that 1 soln per sea state meets criteria
        end

        mult = get_relevant_soln(which_soln,roots,idx_no_sat);   
    end

    assert(all(~isnan(mult),'all'))
end

function mult = get_relevant_soln(which_soln, roots, idx_no_sat)
% pick the specified roots using multidimensional logical indexing

    mult = zeros(size(idx_no_sat));

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

function mult = handle_two_solns(both_ok, which_soln, roots, idx_no_sat,a,b,c)
     
%     % analytical criteria for root between 0 and 1
%     use_1 = (a >= 0 & c <= 0 & (a+b+c) <= 0) | ...
%             (a <= 0 & c <= 0 & (a+b+c) >= 0);
%     use_2 = (a <= 0 & c >= 0 & (a+b+c) <= 0)| ...
%             (a >= 0 & c == 0 & (a+b+c) >= 0);


    [row,col] = find(both_ok);

    % if both are ok, choose the first one arbitrarily for now
    which_soln(row,col,2) = false; 
    mult_1 = get_relevant_soln(which_soln,roots,idx_no_sat);   

    mult = mult_1;
% 
%     % now choose the second one arbitrarily
%     which_soln_2 = which_soln;
%     which_soln_2(row,col,2) = true;
%     which_soln_2(row,col,1) = false;
%     mult_2 = get_relevant_soln(which_soln_2,roots,idx_no_sat);   
% 
%     % compare outliers to figure out if first or second is better
%     [use_1_vals,use_2_vals] = compare_outliers(mult_1,mult_2,both_ok);
%     [use_1_idxs,use_2_idxs] = compare_outliers(double(which_soln(:,:,1)),...
%                                             double(which_soln_2(:,:,1)),both_ok);    
%     
%     use_1 =  (use_1_vals || use_1_idxs) && ~(use_2_vals || use_2_idxs);
%     use_2 = ~(use_1_vals || use_1_idxs) &&  (use_2_vals || use_2_idxs);
% 
%     if use_1
%         mult = mult_1;
%     elseif use_2
%         mult = mult_2;
%     else
%         figure
%         subplot 121
%         contourf(mult_1)
%         subplot 122
%         contourf(mult_2)
% 
%         error(['Failed to figure out which solution to the quadratic ' ...
%             'equation is relevant, try manual inspection.'])
%     end
end

function [use_1,use_2] = compare_outliers(array_1,array_2,relevant_idx)
    window = 6;
    outliers_1 = isoutlier(array_1,'movmedian',window);
    outliers_2 = isoutlier(array_2,'movmedian',window);
    
    outliers_1_relevant = outliers_1(relevant_idx);
    outliers_2_relevant = outliers_2(relevant_idx);

    one_all_outliers = all(outliers_1_relevant);
    two_all_outliers = all(outliers_2_relevant);
    one_all_ok = all(~outliers_1_relevant);
    two_all_ok = all(~outliers_2_relevant);
    
    use_1 = (two_all_outliers && ~one_all_outliers) || (one_all_ok && ~two_all_ok);
    use_2 = (one_all_outliers && ~two_all_outliers) || (two_all_ok && ~one_all_ok);
%     use_1 = sum(outliers_1_relevant,'all') < sum(outliers_2_relevant,'all');
%     use_2 = sum(outliers_1_relevant,'all') > sum(outliers_2_relevant,'all');
end