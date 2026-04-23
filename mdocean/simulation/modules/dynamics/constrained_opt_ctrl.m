function [mag_U,phase_U,...
         real_P,reactive_P,...
         mag_X_u,phase_X_u,...
         mag_X_f,phase_X_f,...
         mag_X_s,phase_X_s,...
         B_p_sat,K_p_sat,...
         qcqp_debug,...
         force_lim_err,...
         amp_lim_err,...
         ctrl_mult_best,...
         phase_ctrl_mult_best,...
         P_unsat] = constrained_opt_ctrl(real_G_u, imag_G_u, w, control_type, ...
                                          control_evaluation_fcn, F_max, P_max, X_max,...
                                          control_solve_type, ctrl_mult_guess, ...
                                          phase_ctrl_mult_guess, T_f_slam, ...
                                          T_s_slam, opt_ctrl_plot_debug_on, ...
                                          B_h_f, B_h_s, B_c, K_h_f, K_h_s, ...
                                          B_drag_f, B_drag_s, m_c, m_f, m_s, ...
                                          multibody, merge_bodies, ...
                                          X_f_guess, X_s_guess, ...
                                          phase_X_f_guess, phase_X_s_guess, ...
                                          H, k_wvn, D_f, D_d)

    % unsaturated optimal control gains
    [B_p,K_p] = controller(real_G_u, imag_G_u, w, control_type);
    stabilize_B = strcmpi(control_type,'reactive') || strcmpi(control_type,'damping');
    stabilize_K = strcmpi(control_type,'reactive');

    % unsaturated response (stabilized) - magnitudes only needed for non-analytical paths
    [mag_U_unsat,P_unsat,...
     mag_X_u_unsat,...
     B_p_stabilized,K_p_stabilized,...
     mag_X_f_unsat,mag_X_s_unsat] = control_evaluation_fcn(K_p,B_p,stabilize_B,stabilize_K);
    
    Z_th = 1 ./ complex(real_G_u, imag_G_u);

    % default qcqp_debug (only populated for 'analytical' path)
    qcqp_debug = struct('centers', [], 'radii', [], 'labels', {{}}, ...
                        'Gamma_opt', NaN, 'alpha', NaN, 'Z_th', NaN, ...
                        'w', NaN, 'n_active_constraints', 0, 'feasible', false);

    if strcmpi(control_type,'none')
        % if no control, no control authority to enforce constraints
        mag_U = mag_U_unsat;
        real_P = P_unsat;
        mag_X_u = mag_X_u_unsat;
        mag_X_f = mag_X_f_unsat;
        mag_X_s = mag_X_s_unsat;
        B_p_sat = B_p_stabilized;
        K_p_sat = K_p_stabilized;

        % phases take a lot longer to compute, so only evaluate them for
        % the no control case - unsat phases aren't needed for the sat case
        [~,real_P,...
         ~,...
         ~,~,...
         ~,~,...
         phase_X_f,phase_X_s,...
         phase_U,reactive_P,...
         phase_X_u] = control_evaluation_fcn(K_p,B_p,stabilize_B,stabilize_K);
        force_lim_err = zeros(size(w));
        amp_lim_err = zeros(size(w));
        % no-control case uses identity multiplier on Z_p (mag=1, phase=0)
        ctrl_mult_best = ones(size(w));
        phase_ctrl_mult_best = zeros(size(w));

    elseif strcmpi(control_solve_type,'analytical')
        % solve constrained optimal control analytically via QCQP circle intersection
        [mag_U_unsat,P_unsat,...
         mag_X_u_unsat,...
         B_p_stabilized,K_p_stabilized,...
         mag_X_f_unsat,mag_X_s_unsat,...
         phase_X_f_unsat,phase_X_s_unsat,...
         phase_U_unsat,~,...
         phase_X_u_unsat] = control_evaluation_fcn(K_p,B_p,true,stabilize_K);
        
        Z_th_complex = 1 ./ (real_G_u + 1i * imag_G_u);

        %Z_f =, Z_s, Z_c to avoid repeat calcs in loop and fewer params passed

        warn_if_infeasible = opt_ctrl_plot_debug_on;
        [B_p_sat,K_p_sat,qcqp_debug] = solve_qcqp_control(Z_th_complex, w, ...
                                mag_X_u_unsat, phase_X_u_unsat, ...
                                mag_X_f_unsat, phase_X_f_unsat, ...
                                mag_X_s_unsat, phase_X_s_unsat, ...
                                mag_U_unsat, phase_U_unsat, ...
                                B_p_stabilized, K_p_stabilized, ...
                                F_max, X_max, ...
                                B_h_f, B_h_s, B_c, K_h_f, K_h_s, ...
                                B_drag_f, B_drag_s, m_c, m_f, m_s, ...
                                control_type, multibody, merge_bodies, ...
                                warn_if_infeasible);

        % evaluate final response with constrained optimal controller
        [mag_U,real_P,...
         mag_X_u,...
         ~,~,...
         mag_X_f,mag_X_s,...
         phase_X_f,phase_X_s,...
         phase_U,reactive_P,...
         phase_X_u] = control_evaluation_fcn(K_p_sat,B_p_sat,true,stabilize_K);

        force_lim_err = zeros(size(w));
        amp_lim_err = zeros(size(w));
        ctrl_mult_best = ones(size(w));
        phase_ctrl_mult_best = zeros(size(w));

    elseif strcmpi(control_solve_type,'solver')
        [mag_U,real_P,...
         mag_X_u,mag_X_f,mag_X_s,...
         B_p_sat,K_p_sat,...
         force_lim_err, ...
         ~,...
         ~,...
         phase_X_f,phase_X_s,...
         amp_lim_err,...
         phase_U,reactive_P,phase_X_u] = response_and_ctrl_err_from_ctrl_guess(...
                                                    ctrl_mult_guess,...
                                                    phase_ctrl_mult_guess,...
                                                    F_max, P_max, mag_U_unsat,...
                                                    X_max, mag_X_f_unsat, ...
                                                    mag_X_s_unsat, ...
                                                    mag_X_u_unsat, ...
                                                    P_unsat, complex(B_p_stabilized, -K_p_stabilized ./ w), Z_th, w,...
                                                    control_evaluation_fcn,...
                                                    H,k_wvn,D_f,D_d,...
                                                    T_f_slam,T_s_slam,...
                                                    control_type);
        ctrl_mult_best = ctrl_mult_guess;
        phase_ctrl_mult_best = phase_ctrl_mult_guess;

    elseif strcmpi(control_solve_type,'brute_force')
        % solve constrained optimal control via brute force
        brute_force_plot_on = opt_ctrl_plot_debug_on;
        [mag_U,phase_U,...
        real_P,reactive_P,...
        mag_X_u,phase_X_u,...
        mag_X_f,phase_X_f,...
        mag_X_s,phase_X_s,...
        B_p_sat,K_p_sat,...
        ctrl_mult_best,phase_ctrl_mult_best] = solve_brute_force_opt_control(Z_th, w, ...
                            mag_X_u_unsat, mag_X_f_unsat, mag_X_s_unsat, ...
                            X_f_guess, phase_X_f_guess, X_s_guess, phase_X_s_guess, ...
                            mag_U_unsat, P_unsat, B_p_stabilized, K_p_stabilized, ...
                             F_max, X_max, P_max, control_type, control_evaluation_fcn, ...
                             H, k_wvn, D_f, D_d, T_f_slam, T_s_slam, brute_force_plot_on && ~merge_bodies);
        force_lim_err = zeros(size(w));
        amp_lim_err = zeros(size(w));
    else
        error('Invalid control_solve_type: %s. Supported values are ''analytical'', ''solver'', or ''brute_force''.', control_solve_type);
    end
end

function [B_p,K_p] = controller(real_G_u, imag_G_u, w, control_type)

    mag_G_u_squared = real_G_u.^2 + imag_G_u.^2;

    % control - set powertrain coefficients
    if strcmpi(control_type,'reactive')
        K_p = -w .* imag_G_u ./ mag_G_u_squared;
        B_p = real_G_u ./ mag_G_u_squared;
    elseif strcmpi(control_type, 'damping')
        B_p = 1 ./ sqrt(mag_G_u_squared);
        K_p = 1e-8; % can't be quite zero because r_k = Inf
    elseif strcmpi(control_type, 'none')
        B_p = zeros(size(w));
        K_p = zeros(size(w));
    end
    B_p(B_p < 0) = 0; % prevent negative controller damping; use indexed assignment to preserve NaN sea states

end

function [B_p_sat,K_p_sat,qcqp_debug] = solve_qcqp_control(Z_th, w, ...
                            mag_X_u_unsat, phase_X_u_unsat, ...
                            mag_X_f_unsat, phase_X_f_unsat, ...
                            mag_X_s_unsat, phase_X_s_unsat, ...
                            mag_U_unsat, phase_U_unsat, ...
                            B_p_stabilized, K_p_stabilized, ...
                            F_max, X_max, ...
                            B_h_f, B_h_s, B_c, K_h_f, K_h_s, ...
                            B_drag_f, B_drag_s, m_c, m_f, m_s, ...
                            control_type, multibody, merge_bodies, ...
                            warn_if_infeasible)
% Solve the constrained optimal control problem analytically using
% the QCQP circle intersection method from the paper.
% Returns the constrained optimal B_p_sat and K_p_sat.

    X_f_max = X_max(1);
    X_s_max = X_max(2);
    X_u_max = X_max(3);
    
    % small tolerance for negligible coefficients (dimensionless and dimensional)
    COEFF_TOL = 1e-10;
    GAMMA_TOL = 1e-12;
    
    B_p_sat = B_p_stabilized;
    K_p_sat = K_p_stabilized;
    qcqp_debug = struct('centers', [], 'radii', [], 'labels', {{}}, ...
                        'Gamma_opt', NaN, 'alpha', NaN, 'Z_th', NaN, ...
                        'w', NaN, 'n_active_constraints', 0, 'feasible', false);
    
    % skip QCQP if no constraints or no PTO
    use_force_sat = isfinite(F_max) && F_max > 0;
    use_amp_sat = any(isfinite(X_max));
    if ~use_force_sat && ~use_amp_sat
        return
    end
    if F_max == 0
        return
    end
    
    % compute complex phasors for the unsaturated (Gamma=0) response
    % and the derivative coefficients c_Y for each quantity Y
    % Y(Gamma) = Y_0 + c_Y * Gamma, where Y_0 is the unsaturated value
    
    for idx = 1:numel(w)
        if isnan(B_h_f(idx))
            continue  % skip NaN sea states
        end
        
        Z_th_i = Z_th(idx);
        w_i = w(idx);
        
        if ~isfinite(Z_th_i) || real(Z_th_i) <= 0
            continue  % skip unstabilizable sea states
        end
        
        % complex phasors at Gamma=0
        X_u_0 = mag_X_u_unsat(idx) * exp(1i * phase_X_u_unsat(idx));
        X_f_0 = mag_X_f_unsat(idx) * exp(1i * phase_X_f_unsat(idx));
        X_s_0 = mag_X_s_unsat(idx) * exp(1i * phase_X_s_unsat(idx));
        U_0   = mag_U_unsat(idx)   * exp(1i * phase_U_unsat(idx));
        
        % I_p = complex current phasor at Gamma=0
        % I = I_p*(1-Gamma), so at Gamma=0: I_0 = I_p
        % velocity = i*w*X_u, so I_p = i*w_i*X_u_0
        I_p = 1i * w_i * X_u_0;
        
        if abs(I_p) < COEFF_TOL
            continue  % no response, nothing to constrain
        end
        
        % derivative coefficients: Y(Gamma) = Y_0 + c_Y * Gamma
        % From the paper: V = I_p*Z_th^**(1+Gamma), I = I_p*(1-Gamma)
        c_V = I_p * conj(Z_th_i);   % force derivative: dV/dGamma
        c_I = -I_p;                   % velocity derivative: dI/dGamma
        c_Xu = -X_u_0;               % PTO displacement: X_u = I/(iw)
        
        % float and spar displacement derivatives depend on topology
        if ~multibody || merge_bodies
            % single body: X_f = X_u, X_s = 0
            c_Xf = c_Xu;
            c_Xs = 0;
        else
            % multibody: compute transfer functions from impedance matrices
            % All arrays indexed by (idx) to extract the scalar for this sea state
            B_f_total = B_h_f(idx) + B_drag_f(idx);
            B_s_total = B_h_s(idx) + B_drag_s(idx);
            Z_f_i = B_f_total + 1i*(m_f(idx)*w_i - K_h_f(idx)/w_i);
            Z_s_i = B_s_total + 1i*(m_s(idx)*w_i - K_h_s(idx)/w_i);
            Z_c_i = B_c(idx) + 1i*(m_c(idx)*w_i);
            det_Z_i = Z_f_i .* Z_s_i - Z_c_i.^2;
            
            % transfer from PTO force to body displacement
            % v_f = v_f_forced + (Z_s+Z_c)/det_Z * F_pto
            % X_f = X_f_forced + (Z_s+Z_c)/(det_Z*iw) * F_pto
            % F_pto = V = I_p * Z_th^* * (1+Gamma)
            % c_Xf = (Z_s+Z_c)/(det_Z*iw) * I_p * Z_th^*
            T_f = (Z_s_i + Z_c_i) ./ (det_Z_i .* 1i .* w_i);
            T_s = -(Z_f_i + Z_c_i) ./ (det_Z_i .* 1i .* w_i);
            c_Xf = T_f .* I_p .* conj(Z_th_i);
            c_Xs = T_s .* I_p .* conj(Z_th_i);
        end
        
        % build circle constraints
        % For constraint |Y| <= Y_max where Y = Y_0 + c_Y * Gamma:
        % circle center = -Y_0/c_Y, radius = Y_max/|c_Y|, S=+1 (inside)
        centers = [];
        radii = [];
        current_labels = {};
        
        % force constraint: |V| <= F_max
        if use_force_sat && abs(c_V) > COEFF_TOL
            center_F = -U_0 / c_V;
            radius_F = F_max / abs(c_V);
            centers = [centers; real(center_F), imag(center_F)];
            radii = [radii; radius_F];
            current_labels{end+1} = 'Force limit';
        end
        
        % float amplitude constraint: |X_f| <= X_f_max
        if isfinite(X_f_max) && abs(c_Xf) > COEFF_TOL
            center_Xf = -X_f_0 / c_Xf;
            radius_Xf = X_f_max / abs(c_Xf);
            centers = [centers; real(center_Xf), imag(center_Xf)];
            radii = [radii; radius_Xf];
            current_labels{end+1} = 'Float amplitude';
        end
        
        % spar amplitude constraint: |X_s| <= X_s_max
        if isfinite(X_s_max) && abs(c_Xs) > COEFF_TOL && multibody && ~merge_bodies
            center_Xs = -X_s_0 / c_Xs;
            radius_Xs = X_s_max / abs(c_Xs);
            centers = [centers; real(center_Xs), imag(center_Xs)];
            radii = [radii; radius_Xs];
            current_labels{end+1} = 'Spar amplitude';
        end
        
        % PTO amplitude constraint: |X_u| <= X_u_max
        if isfinite(X_u_max) && abs(c_Xu) > COEFF_TOL
            center_Xu = -X_u_0 / c_Xu;
            radius_Xu = X_u_max / abs(c_Xu);
            
            centers = [centers; real(center_Xu), imag(center_Xu)];
            radii = [radii; radius_Xu];
            current_labels{end+1} = 'PTO amplitude';
        end
        
        % positive power constraint: |Gamma - i*alpha|^2 <= 1 + alpha^2
        % where alpha = Im(Z_th)/Re(Z_th)
        alpha = imag(Z_th_i) / real(Z_th_i);
        center_P = [0, alpha];
        radius_P = sqrt(1 + alpha^2);
        % always include positive power constraint (it's cheap and prevents P<0)
        centers = [centers; center_P];
        radii = [radii; radius_P];
        current_labels{end+1} = 'Positive power';
        
        % damping control: add Q=0 constraint circle
        % Gamma must lie ON the Q=0 circle (equality constraint)
        % Q=0 circle: center = -i*R/X, radius = |Z_th|/|X|
        % where R = Re(Z_th), X = Im(Z_th)
        p_star = [];
        if strcmpi(control_type, 'damping') && abs(imag(Z_th_i)) > COEFF_TOL
            R_th = real(Z_th_i);
            X_th = imag(Z_th_i);
            center_Q0 = [0, -R_th/X_th];
            radius_Q0 = abs(Z_th_i) / abs(X_th);
            % for damping, find optimal on Q=0 circle subject to other constraints
            Gamma_opt = solve_damping_qcqp(center_Q0, radius_Q0, centers, radii);
        else
            % reactive control or no reactance: use standard circle intersection
            if isempty(centers)
                Gamma_opt = 0;  % no constraints active
            else
                [p_star, ~, ~] = circle_intersect_optim(centers, radii, warn_if_infeasible);
                if isempty(p_star)
                    % Constraint circles are mutually infeasible.
                    % This typically happens when force (circle at Gamma=-1)
                    % and amplitude (circle at Gamma=+1) limits are both
                    % violated and r_F + r_Xu < 2 (they can't intersect).
                    % Fall back: try each individual constraint circle with
                    % the positive-power circle to find best feasible Gamma.
                    % Power circle is always the last entry in centers/radii.
                    Gamma_opt = 1;  % default: shutdown
                    best_norm = Inf;
                    n_c = size(centers, 1);
                    c_pow = centers(n_c, :);
                    r_pow = radii(n_c);
                    for k = 1:n_c - 1
                        % closest boundary point on circle k to origin
                        dist_to_center_k = norm(centers(k,:));
                        if dist_to_center_k < COEFF_TOL
                            boundary_point_k = [radii(k), 0];
                        else
                            boundary_point_k = centers(k,:) - radii(k) * centers(k,:)/dist_to_center_k;
                        end
                        % accept only if inside the positive-power circle
                        if norm(boundary_point_k - c_pow) <= r_pow + 1e-4 && norm(boundary_point_k) < best_norm
                            best_norm = norm(boundary_point_k);
                            Gamma_opt = boundary_point_k(1) + 1i*boundary_point_k(2);
                        end
                    end
                else
                    Gamma_opt = p_star(1) + 1i*p_star(2);
                end
            end
        end
        
        % update debug struct if this is the most constrained feasible sea state
        n_non_power = size(centers, 1) - 1;  % last circle is always the power circle
        if n_non_power > qcqp_debug.n_active_constraints && ~isempty(centers) && isfinite(Gamma_opt)
            qcqp_debug.centers = centers;
            qcqp_debug.radii = radii;
            qcqp_debug.labels = current_labels;
            qcqp_debug.Gamma_opt = Gamma_opt;
            qcqp_debug.alpha = imag(Z_th_i) / real(Z_th_i);
            qcqp_debug.Z_th = Z_th_i;
            qcqp_debug.w = w_i;
            qcqp_debug.n_active_constraints = n_non_power;
            qcqp_debug.feasible = ~isempty(p_star);
        end

        % convert Gamma to Z_l, then to B_p and K_p
        if abs(Gamma_opt) < GAMMA_TOL
            % no constraint active, keep unsaturated gains
            continue
        end
        
        Z_l = conj(Z_th_i) * (1 + Gamma_opt) / (Gamma_opt - 1);
        B_p_sat(idx) = real(Z_l);
        K_p_sat(idx) = -w_i * imag(Z_l);
        
        % ensure B_p >= 0 (stabilize)
        if B_p_sat(idx) < 0
            B_p_sat(idx) = 0;
        end
    end
end

function Gamma_opt = solve_damping_qcqp(center_Q0, radius_Q0, centers, radii)
% For damping control, find the point on the Q=0 circle closest to the 
% origin that satisfies all other constraint circles.
% The Q=0 circle parametrically: Gamma(t) = center + radius * exp(i*t)

    % closest point on Q=0 circle to origin (unconstrained damping optimal)
    c_Q0 = center_Q0(1) + 1i*center_Q0(2);
    dist_to_center = abs(c_Q0);
    if dist_to_center < 1e-10
        % center at origin, any point on circle is equally close
        Gamma_closest = radius_Q0;
    else
        Gamma_closest = c_Q0 - radius_Q0 * c_Q0 / dist_to_center;
    end
    
    % check if the closest point is feasible
    p_closest = [real(Gamma_closest), imag(Gamma_closest)];
    if isempty(centers) || check_inside_circles(p_closest, centers, radii)
        Gamma_opt = Gamma_closest;
        return
    end
    
    % find intersections of Q=0 circle with each constraint circle
    candidates = [];
    for i = 1:length(radii)
        [xs, ys] = circle_circle_intersect_inline(...
            center_Q0(1), center_Q0(2), radius_Q0, ...
            centers(i,1), centers(i,2), radii(i));
        if ~isnan(xs(1))
            candidates = [candidates; xs(1), ys(1); xs(2), ys(2)]; %#ok<AGROW>
        end
    end
    
    % filter feasible candidates
    feasible = [];
    for i = 1:size(candidates,1)
        p = candidates(i,:);
        % check on Q=0 circle and inside all other constraint circles
        on_Q0 = abs(norm(p - center_Q0) - radius_Q0) < 1e-4;
        if on_Q0 && check_inside_circles(p, centers, radii)
            feasible = [feasible; p]; %#ok<AGROW>
        end
    end
    
    if isempty(feasible)
        % No point on the Q=0 circle satisfies all constraints simultaneously.
        % Fall back to the unconstrained damping-optimal point (Gamma_closest),
        % which always satisfies the positive-power constraint (the Q=0 and
        % power circles always intersect when Im(Z_th)≠0).
        % This may violate force/amplitude limits but avoids complete shutdown.
        Gamma_opt = Gamma_closest;
    else
        % choose closest to origin
        dists = vecnorm(feasible, 2, 2);
        [~, best] = min(dists);
        Gamma_opt = feasible(best,1) + 1i*feasible(best,2);
    end
end

function inside = check_inside_circles(p, centers, radii)
% Check if point p=[x,y] is inside all circles
    dists = vecnorm(p - centers, 2, 2);
    inside = all(dists <= radii + 1e-4);
end

function [xs, ys] = circle_circle_intersect_inline(x1,y1,r1,x2,y2,r2)
% Inline circle-circle intersection (same as in circle_intersect_optim)
    d = sqrt((x2-x1)^2 + (y2-y1)^2);
    if d > r1+r2 || d < abs(r1-r2) || d == 0
        xs = [NaN, NaN];
        ys = [NaN, NaN];
        return
    end
    a = (r1^2 - r2^2 + d^2) / (2*d);
    h = sqrt(max(r1^2 - a^2, 0));
    mx = x1 + a*(x2-x1)/d;
    my = y1 + a*(y2-y1)/d;
    dx = h*(y2-y1)/d;
    dy = h*(x2-x1)/d;
    xs = [mx+dx, mx-dx];
    ys = [my-dy, my+dy];
end

function [opt_mag_U,opt_phase_U,...
    opt_real_P,opt_reactive_P,...
    opt_mag_X_u,opt_phase_X_u,...
    opt_mag_X_f,opt_phase_X_f,...
    opt_mag_X_s,opt_phase_X_s,...
    opt_B_p_sat,opt_K_p_sat,...
    opt_ctrl_mult,opt_phase_ctrl_mult] = solve_brute_force_opt_control(Z_th, w, ...
                            mag_X_u_unsat, mag_X_f_unsat, mag_X_s_unsat, ...
                            X_f_guess, phase_X_f_guess, X_s_guess, phase_X_s_guess, ...
                            mag_U_unsat, P_unsat, B_p_unsat, K_p_unsat, ...
                            F_max, X_max, P_max, control_type, ...
                            control_evaluation_fcn, H, k_wvn, D_f, D_d,...
                            T_f_slam,T_s_slam, brute_force_plot_on)
    
    Z_p_unsat = complex(B_p_unsat, - K_p_unsat ./ w);

    size_opt_ctrl_mesh = 11;  
    if rem(size_opt_ctrl_mesh,2)==0 % input is even
        size_mult_mesh  = size_opt_ctrl_mesh + 1; % this MUST be odd so ctrl_mult_guess includes 1!
        size_phase_mesh = size_opt_ctrl_mesh;    % this should ideally be even for symmetrical plots
    else
        size_mult_mesh  = size_opt_ctrl_mesh;
        size_phase_mesh = size_opt_ctrl_mesh + 1;
    end
    ctrl_mult_guess = logspace(-1,1,size_mult_mesh);
    if strcmp(control_type,'reactive')
        phase_ctrl_mult_guess = linspace(0,2*pi,size_phase_mesh+1);
        phase_ctrl_mult_guess = phase_ctrl_mult_guess(1:end-1); % 2pi=0 so don't repeat it
    else
        phase_ctrl_mult_guess = 0;
    end
    [MAG_GUESS,PHASE_GUESS] = meshgrid(ctrl_mult_guess,phase_ctrl_mult_guess);

    mag_U = nan([size(w), numel(MAG_GUESS)]);
    real_P = nan([size(w), numel(MAG_GUESS)]);
    mag_X_u = nan([size(w), numel(MAG_GUESS)]);
    mag_X_f = nan([size(w), numel(MAG_GUESS)]);
    mag_X_s = nan([size(w), numel(MAG_GUESS)]);
    B_p_sat = nan([size(w), numel(MAG_GUESS)]);
    K_p_sat = nan([size(w), numel(MAG_GUESS)]);
    constraint_err = nan([size(w), numel(MAG_GUESS)]);
    mag_ctrl_mult_stabilized = nan([size(w), numel(MAG_GUESS)]);
    phase_ctrl_mult_stabilized = nan([size(w), numel(MAG_GUESS)]);

    for i = 1:numel(MAG_GUESS)
        [mag_U(:,:,i),real_P_out(:,:,i),...
        mag_X_u(:,:,i),mag_X_f(:,:,i),mag_X_s(:,:,i),...
        B_p_sat(:,:,i),K_p_sat(:,:,i),...
        constraint_err(:,:,i),...
        mag_ctrl_mult_stabilized(:,:,i),...
        phase_ctrl_mult_stabilized(:,:,i)] = response_and_ctrl_err_from_ctrl_guess(MAG_GUESS(i)*ones(size(w)),...
                                                    PHASE_GUESS(i)*ones(size(w)),...
                                                    F_max, P_max, mag_U_unsat,...
                                                    X_max, mag_X_f_unsat, ...
                                                    mag_X_s_unsat, ...
                                                    mag_X_u_unsat, ...
                                                    P_unsat, Z_p_unsat, Z_th, w,...
                                                    control_evaluation_fcn,...
                                                    H,k_wvn,D_f,D_d,...
                                                    T_f_slam,T_s_slam,...
                                                    control_type,...
                                                    X_f_guess, phase_X_f_guess,...
                                                    X_s_guess, phase_X_s_guess);
        if strcmpi(control_type,'reactive')
            P_matched = P_unsat;
            z = mag_ctrl_mult_stabilized(:,:,i) .* exp( 1i * phase_ctrl_mult_stabilized(:,:,i) );
        else
            error('brute force opt ctrl not yet supported for non reactive control')
        end
        mag_gamma = abs((z-1)./(z+1));
        real_P(:,:,i) = P_matched .* (1 - mag_gamma.^2);
    end

    % save unmodified power for polar plot before infeasible controllers are set to -1
    real_P_unmodified = real_P;

    % choose best controller for each sea state
    each_controller_feasible = constraint_err <= 0;
    ctrl_dim = 3;
    sea_state_dim =[1 2];
    each_sea_state_feasible = any(each_controller_feasible,ctrl_dim);
    real_P(each_sea_state_feasible & ~each_controller_feasible) = -1; % negative power for infeasible controllers in sea states where a feasbile controller exists, so they aren't chosen
    [opt_real_P,idx_opt] = max(real_P,[],ctrl_dim);

    sea_state_infeasible = ~each_sea_state_feasible & ~isnan(w);
    if any(sea_state_infeasible(:))
        
        constraint_err_pos_pwr = constraint_err;
        constraint_err_pos_pwr(real_P_unmodified<0) = Inf;
        
        [min_err,idx_least_err] = min(constraint_err_pos_pwr,[],ctrl_dim);
        idx_opt(sea_state_infeasible) = idx_least_err(sea_state_infeasible);
        % Use sub2ind for correct 3D indexing: idx_opt gives ctrl-dim index for
        % each (Hs,T) sea state, so plain linear indexing real_P_unmodified(idx_opt)
        % would silently access the wrong elements of the [Nh x NT x N_ctrl] array.
        sz = size(real_P_unmodified);
        [r,c] = ndgrid(1:sz(1), 1:sz(2));
        lin_idx = sub2ind(sz, r, c, idx_opt);
        opt_real_P = real_P_unmodified(lin_idx);

        debug_print = brute_force_plot_on;
        if debug_print
            warning(['no feasible controller found for %d/%d sea states, setting signals to '...
                'least error solution. You may want to adjust the ctrl_mult_guess to be wider'], ...
                sum(sea_state_infeasible(:)), sum(~isnan(w(:))))
            % diagnostic: identify which constraint(s) are not satisfied
            infeas_mask = ~each_sea_state_feasible;
            fprintf('  Brute force diagnostics for infeasible sea states:\n');
            fprintf('    min constraint_err across controllers: min=%.4g, max=%.4g\n', ...
                    min(min_err(infeas_mask)), max(min_err(infeas_mask)));
    
            % diagnostic: identify the number of sea states for which each
            % controller was able to satisfy all constraints
            num_ss_ctrl_feasible = reshape(squeeze(sum(each_controller_feasible,sea_state_dim)),size(MAG_GUESS));
            mag_phs_num_feasible = [MAG_GUESS(:) PHASE_GUESS(:) num_ss_ctrl_feasible(:)];
            fprintf('   Controllers and the number of sea states they are feasible: \n[mag,phase,number]=\n[%s]\n', ...
                formattedDisplayText(mag_phs_num_feasible))
    
            % evaluate per-constraint errors at the best guess for each infeasible sea state
            best_B_p = nan(size(w));
            best_mag_U_val = nan(size(w));
            best_mag_X_f_val = nan(size(w));
            best_mag_X_s_val = nan(size(w));
            best_mag_X_u_val = nan(size(w));
            best_real_P_val = nan(size(w));
            for idx_flat = find(infeas_mask(:))'
                linear_idx = idx_flat + (idx_least_err(idx_flat)-1)*numel(w);
                best_B_p(idx_flat) = B_p_sat(linear_idx);
                best_mag_U_val(idx_flat) = mag_U(linear_idx);
                best_mag_X_f_val(idx_flat) = mag_X_f(linear_idx);
                best_mag_X_s_val(idx_flat) = mag_X_s(linear_idx);
                best_mag_X_u_val(idx_flat) = mag_X_u(linear_idx);
                best_real_P_val(idx_flat) = real_P(linear_idx);
            end
            n_infeas = sum(infeas_mask(:));
            fprintf('    B_p<0 (neg damping):   %d/%d infeasible sea states\n', sum(best_B_p(infeas_mask) < 0), n_infeas);
            fprintf('    force exceeded F_max:  %d/%d (F_max=%.4g)\n', sum(best_mag_U_val(infeas_mask) > 4/pi*F_max), n_infeas, F_max);
            fprintf('    power exceeded P_max:  %d/%d (P_max=%.4g)\n', sum(best_real_P_val(infeas_mask) > P_max), n_infeas, P_max);
            fprintf('    float amp exceeded:    %d/%d (X_max(1)=%.4g)\n', sum(best_mag_X_f_val(infeas_mask) > X_max(1)), n_infeas, X_max(1));
            fprintf('    spar amp exceeded:     %d/%d (X_max(2)=%.4g)\n', sum(best_mag_X_s_val(infeas_mask) > X_max(2)), n_infeas, X_max(2));
            fprintf('    rel amp exceeded:      %d/%d (X_max(3)=%.4g)\n', sum(best_mag_X_u_val(infeas_mask) > X_max(3)), n_infeas, X_max(3));
        end
    end

    % polar plot of controller search grid (only in operational, non-mergedbodies case)
    if brute_force_plot_on
        make_ctrl_polar_plot(MAG_GUESS, PHASE_GUESS, real_P_unmodified, ...
                                constraint_err, idx_opt, ...
                                mag_ctrl_mult_stabilized, phase_ctrl_mult_stabilized);
    end

    % Re-evaluate the optimal controller once to obtain all phase outputs that were
    % skipped during the brute-force loop (phase_U, reactive_P, phase_X_u, phase_X_f,
    % phase_X_s are not needed to find the optimum and are only computed here).
    [opt_mag_U, ~, opt_mag_X_u, opt_mag_X_f, opt_mag_X_s, ...
     opt_B_p_sat, opt_K_p_sat, ~, opt_ctrl_mult, opt_phase_ctrl_mult, ...
     opt_phase_X_f, opt_phase_X_s, ...
     ~, opt_phase_U, opt_reactive_P, opt_phase_X_u] = ...
        response_and_ctrl_err_from_ctrl_guess(MAG_GUESS(idx_opt), PHASE_GUESS(idx_opt), ...
                                              F_max, P_max, mag_U_unsat, X_max, ...
                                              mag_X_f_unsat, mag_X_s_unsat, mag_X_u_unsat, ...
                                              P_unsat, Z_p_unsat, Z_th, w, ...
                                              control_evaluation_fcn, H, k_wvn, D_f, D_d, ...
                                              T_f_slam, T_s_slam, control_type);

    % Warn if the final optimal response contains NaN for valid sea states.
    % This path is not covered by the warning in dynamics_error_wrapper because
    % that function is only called during the iterative drag-convergence loop.
    if any(isnan(opt_mag_X_f(~isnan(w))),'all')
        warning('NaN in opt_mag_X_f after brute-force re-evaluation for non-NaN sea states')
    end
    if any(isnan(opt_real_P(~isnan(w))),'all')
        warning('NaN in opt_real_P after brute-force re-evaluation for non-NaN sea states')
    end
end

function [mag_U,real_P,...
        mag_X_u,mag_X_f,mag_X_s,...
        B_stabilized,K_stabilized,...
        constr_viol_err,...
        mag_ctrl_mult_stabilized,...
        phase_ctrl_mult_stabilized,...
        phase_X_f,phase_X_s,...
        optimality_err,...
        phase_U,reactive_P,phase_X_u] = response_and_ctrl_err_from_ctrl_guess(...
                                                    ctrl_mult_guess,...
                                                    phase_ctrl_mult_guess,...
                                                    F_max, P_max, mag_U_unsat,...
                                                    X_max, mag_X_f_unsat, ...
                                                    mag_X_s_unsat, ...
                                                    mag_X_u_unsat, ...
                                                    real_P_unsat,Z_p, Z_th, w,...
                                                    control_evaluation_fcn,...
                                                    H,k,D_f,D_d,T_f,T_s,...
                                                    control_type,...
                                                    mag_X_f_guess, phase_X_f_guess,...
                                                    mag_X_s_guess, phase_X_s_guess)

    % saturated control gains from guess
    mult = ctrl_mult_guess .* exp(1i*phase_ctrl_mult_guess);
    Z_p_sat = Z_p .* mult;
    B_p_sat =       real(Z_p_sat);
    K_p_sat = -w .* imag(Z_p_sat);

    % saturated response (with B stabilization, and K stabilization for reactive control)
    stabilize_B = true;
    stabilize_K = strcmpi(control_type,'reactive');

    % Determine whether pre-computed X_f/X_s guesses were supplied (fast path).
    % When nargin > 21 the caller passes guessed mag/phase for X_f and X_s so
    % that the brute-force loop only needs to evaluate power, mag_U, and mag_X_u
    % from the controller (outputs 1-5 of evaluate_any_controller after the
    % B_l/K_l reorder).  The guesses avoid calling multibody_response for the
    % X_f/X_s rows and their phases in every loop iteration.
    use_X_guess = nargin > 21;

    if use_X_guess
        % Fast path: only evaluate power, mag_U, mag_X_u, B_stabilized, K_stabilized.
        % Outputs 1-5 of evaluate_any_controller after the B_l/K_l reorder are
        % [mag_U, real_P, mag_X_u, B_l, K_l], so nargout=5 skips X_f/X_s.
        [mag_U,real_P,...
         mag_X_u,...
         B_stabilized,K_stabilized] = control_evaluation_fcn(K_p_sat,B_p_sat,stabilize_B,stabilize_K);
        mag_X_f  = mag_X_f_guess;
        phase_X_f = phase_X_f_guess;
        mag_X_s  = mag_X_s_guess;
        phase_X_s = phase_X_s_guess;
    elseif nargout > 13
        % caller needs phase_U, reactive_P, phase_X_u — request all outputs
        [mag_U,real_P,...
         mag_X_u,...
         B_stabilized,K_stabilized,...
         mag_X_f,mag_X_s,...
         phase_X_f,phase_X_s,...
         phase_U,reactive_P,...
         phase_X_u] = control_evaluation_fcn(K_p_sat,B_p_sat,stabilize_B,stabilize_K);
    else
        % caller only needs magnitudes and phase_X_f/phase_X_s (for slamming constraint)
        [mag_U,real_P,...
         mag_X_u,...
         B_stabilized,K_stabilized,...
         mag_X_f,mag_X_s,...
         phase_X_f,phase_X_s] = control_evaluation_fcn(K_p_sat,B_p_sat,stabilize_B,stabilize_K);
    end
    Z_p_stabilized = complex(B_stabilized, -K_stabilized ./ w);

    stabilization_mult = Z_p_stabilized ./ Z_p_sat;
    ctrl_mult_stabilized = mult .* stabilization_mult;
    mag_ctrl_mult_stabilized = abs(ctrl_mult_stabilized);
    phase_ctrl_mult_stabilized = my_angle(ctrl_mult_stabilized);

    % control error: uses saturated and unsaturated response
    if nargout > 12
        [constr_viol_err,optimality_err] = control_errors_from_sat_results(ctrl_mult_guess,phase_ctrl_mult_guess,...
                                                    F_max,P_max,mag_U_unsat,...
                                                    X_max(1), mag_X_f_unsat, mag_X_f, phase_X_f,...
                                                    X_max(2), mag_X_s_unsat, mag_X_s, phase_X_s,...
                                                    X_max(3), mag_X_u_unsat,...
                                                    real_P_unsat, Z_th, ...
                                                    H, k, D_f, D_d, T_f, T_s,...
                                                    mag_U, mag_X_u, B_p_sat, real_P);
    else
        constr_viol_err = control_errors_from_sat_results(ctrl_mult_guess,phase_ctrl_mult_guess,...
                                                    F_max,P_max,mag_U_unsat,...
                                                    X_max(1), mag_X_f_unsat, mag_X_f, phase_X_f,...
                                                    X_max(2), mag_X_s_unsat, mag_X_s, phase_X_s,...
                                                    X_max(3), mag_X_u_unsat,...
                                                    real_P_unsat, Z_th, ...
                                                    H, k, D_f, D_d, T_f, T_s,...
                                                    mag_U, mag_X_u, B_p_sat, real_P);
    end
end


function [constr_viol_err,optimality_err] = control_errors_from_sat_results(ctrl_mult_mag,...
                                                   ctrl_mult_phase, F_max, P_max, mag_U_unsat,...
                                                   X_f_upper_limit_static, mag_X_f_unsat, mag_X_f, phase_X_f,...
                                                   X_s_upper_limit_static, mag_X_s_unsat, mag_X_s, phase_X_s,...
                                                   X_u_upper_limit_static, mag_X_u_unsat,...
                                                   P_unsat, Z_th,...
                                                   H, k, D_f, D_d, T_f, T_s,...
                                                   mag_U, mag_X_u, B_p_sat, P_sat)

    use_amp_sat = ~(isinf(X_f_upper_limit_static) && isinf(X_s_upper_limit_static) && isinf(X_u_upper_limit_static));

    if use_amp_sat
        [X_f_lower_limit_dynamic,...
         X_f_upper_limit_dynamic,...
         idx_f_imag] = get_slamming_min_max(H/2, k, D_f, phase_X_f, T_f);
        X_f_lower_limit_dynamic(idx_f_imag) = -abs(imag(X_f_lower_limit_dynamic(idx_f_imag)));
        X_f_upper_limit_dynamic(idx_f_imag) = -abs(imag(X_f_upper_limit_dynamic(idx_f_imag)));

        [X_s_lower_limit_dynamic,...
         X_s_upper_limit_dynamic,...
         idx_s_imag] = get_slamming_min_max(H/2, k, D_d, phase_X_s, T_s);
        X_s_upper_limit_dynamic(idx_s_imag) = -abs(imag(X_s_upper_limit_dynamic(idx_s_imag)));
        X_s_lower_limit_dynamic(idx_s_imag) = -abs(imag(X_s_lower_limit_dynamic(idx_s_imag)));

        X_f_upper_limit = min(X_f_upper_limit_dynamic, X_f_upper_limit_static);
        X_s_upper_limit = min(X_s_upper_limit_dynamic, X_s_upper_limit_static);
        X_u_upper_limit = X_u_upper_limit_static;

        X_f_lower_limit = X_f_lower_limit_dynamic;
        X_s_lower_limit = X_s_lower_limit_dynamic;
    else
        X_f_upper_limit = Inf(size(mag_X_f));
        X_s_upper_limit = Inf(size(mag_X_s));
        X_u_upper_limit = Inf(size(mag_X_u));
        X_f_lower_limit = -Inf(size(mag_X_f));
        X_s_lower_limit = -Inf(size(mag_X_s));
    end

    % what the force, amp, and power would be if the force-sat, amp-sat,
    % and power-sat solutions applied - notebook p106 2/2/25
    force_saturated_U = min(4/pi * F_max,      mag_U_unsat);

    amp_up_saturated_X_f = min(X_f_upper_limit, mag_X_f_unsat);
    amp_dn_saturated_X_f = max(X_f_lower_limit, mag_X_f_unsat);

    amp_up_saturated_X_s = min(X_s_upper_limit, mag_X_s_unsat);
    amp_dn_saturated_X_s = max(X_s_lower_limit, mag_X_s_unsat);

    amp_saturated_X_u = min(X_u_upper_limit, mag_X_u_unsat);

    P_min = 0;
    power_up_saturated_P = min(P_max, P_unsat);
    power_dn_saturated_P = max(P_min, P_unsat);

    % indices where sat and unsat solutions violate which constraint
    idx_force_viol_unsat    = mag_U_unsat   > force_saturated_U;
    idx_amp_f_viol_unsat_up = mag_X_f_unsat > amp_up_saturated_X_f;
    idx_amp_f_viol_unsat_dn = mag_X_f_unsat < amp_dn_saturated_X_f;
    idx_amp_s_viol_unsat_up = mag_X_s_unsat > amp_up_saturated_X_s;
    idx_amp_s_viol_unsat_dn = mag_X_s_unsat < amp_dn_saturated_X_s;
    idx_amp_u_viol_unsat    = mag_X_u_unsat > amp_saturated_X_u;
    idx_power_viol_unsat_up = P_unsat       > P_max;
    idx_power_viol_unsat_dn = P_unsat       < P_min;
    
    num_constr_viol_unsat = idx_force_viol_unsat    + idx_amp_f_viol_unsat_up ...
                          + idx_amp_f_viol_unsat_dn + idx_amp_s_viol_unsat_up ...
                          + idx_amp_s_viol_unsat_dn + idx_amp_u_viol_unsat ...
                          + idx_power_viol_unsat_up + idx_power_viol_unsat_dn;
    all_ok_unsat       = num_constr_viol_unsat==0 & ~isnan(P_unsat);
    only_force_viol_unsat    = idx_force_viol_unsat    & num_constr_viol_unsat==1;
    only_amp_f_viol_unsat_up = idx_amp_f_viol_unsat_up & num_constr_viol_unsat==1;
    only_amp_f_viol_unsat_dn = idx_amp_f_viol_unsat_dn & num_constr_viol_unsat==1;
    only_amp_s_viol_unsat_up = idx_amp_s_viol_unsat_up & num_constr_viol_unsat==1;
    only_amp_s_viol_unsat_dn = idx_amp_s_viol_unsat_dn & num_constr_viol_unsat==1;
    only_amp_u_viol_unsat    = idx_amp_u_viol_unsat    & num_constr_viol_unsat==1;
    only_power_viol_unsat_up = idx_power_viol_unsat_up & num_constr_viol_unsat==1;
    only_power_viol_unsat_dn = idx_power_viol_unsat_dn & num_constr_viol_unsat==1;
    mult_const_viol_unsat = num_constr_viol_unsat >= 2;

    % indices where each solution applies. commented forumlas attempt to 
    % incorporate both but should not actually be based purely on whether
    % the saturated and/or unsaturated solution violates, it should also
    % depend on alpha and the limits (see notebook p142-144 9/15/25).
    idx_force_sat_applies    = all_ok_unsat | only_force_viol_unsat;
    idx_amp_f_sat_applies_up = all_ok_unsat | only_amp_f_viol_unsat_up;
    idx_amp_f_sat_applies_dn = all_ok_unsat | only_amp_f_viol_unsat_dn;
    idx_amp_s_sat_applies_up = all_ok_unsat | only_amp_s_viol_unsat_up;
    idx_amp_s_sat_applies_dn = all_ok_unsat | only_amp_s_viol_unsat_dn;
    idx_amp_u_sat_applies    = all_ok_unsat | only_amp_u_viol_unsat;
    idx_power_sat_applies_up = all_ok_unsat | only_power_viol_unsat_up;
    idx_power_sat_applies_dn = all_ok_unsat | only_power_viol_unsat_dn;

    % errors from each solution: >0 when constr is violated, and <0
    % when constr is not tight. 
    if isinf(F_max)
        F_err_from_force_sat = zeros(size(mag_X_f));
    else
        F_err_from_force_sat = mag_U   ./ force_saturated_U - 1;
        idx_U_0 = force_saturated_U == 0;
        F_err_from_force_sat(idx_U_0) = mag_U(idx_U_0) / 1e6;
        F_err_from_force_sat(F_err_from_force_sat < -1) = -.99;
    end
    X_err_from_amp_f_sat_up = mag_X_f ./ amp_up_saturated_X_f - 1;
    idx_amp_f_0 = amp_up_saturated_X_f == 0;
    X_err_from_amp_f_sat_up(idx_amp_f_0) = mag_X_f(idx_amp_f_0);
    X_err_from_amp_f_sat_up(X_err_from_amp_f_sat_up < -1) = -.99;

    X_err_from_amp_f_sat_dn = amp_dn_saturated_X_f ./ mag_X_f  - 1;
    idx_X_f_0 = mag_X_f == 0;
    X_err_from_amp_f_sat_dn(idx_X_f_0) = -amp_dn_saturated_X_f(idx_X_f_0);
    X_err_from_amp_f_sat_dn(X_err_from_amp_f_sat_dn < -1) = -.99;

    X_err_from_amp_s_sat_up = mag_X_s ./ amp_up_saturated_X_s - 1;
    idx_amp_s_0 = amp_up_saturated_X_s == 0;
    X_err_from_amp_s_sat_up(idx_amp_s_0) = mag_X_s(idx_amp_s_0);
    X_err_from_amp_s_sat_up(X_err_from_amp_s_sat_up < -1) = -.99;

    X_err_from_amp_s_sat_dn = amp_dn_saturated_X_s ./ mag_X_s  - 1;
    idx_X_s_0 = mag_X_s == 0;
    X_err_from_amp_s_sat_dn(idx_X_s_0) = -amp_dn_saturated_X_s(idx_X_s_0);
    X_err_from_amp_s_sat_dn(X_err_from_amp_s_sat_dn < -1) = -.99;

    X_err_from_amp_u_sat = mag_X_u ./ amp_saturated_X_u - 1;
    idx_amp_u_0 = amp_saturated_X_u == 0;
    X_err_from_amp_u_sat(idx_amp_u_0) = mag_X_u(idx_amp_u_0);
    X_err_from_amp_u_sat(X_err_from_amp_u_sat < -1) = -.99;

    if isinf(P_max)
        P_err_from_power_sat_up = zeros(size(mag_X_f));
    else
        P_err_from_power_sat_up = P_sat ./ power_up_saturated_P - 1;
        idx_P_up_0 = power_up_saturated_P == 0;
        P_err_from_power_sat_up(idx_P_up_0) = P_sat(idx_P_up_0) / 1e6;
        P_err_from_power_sat_up(P_err_from_power_sat_up < -1) = -.99;
    end

    P_err_from_power_sat_dn = power_dn_saturated_P ./ P_sat - 1;
    idx_P_dn_0 = P_sat == 0;
    P_err_from_power_sat_dn(idx_P_dn_0) = -power_dn_saturated_P(idx_P_dn_0) / 1e6;
    P_err_from_power_sat_dn(P_err_from_power_sat_dn < -1) = -.99;

    if ~use_amp_sat
        X_err_from_amp_f_sat_up = zeros(size(mag_X_f));
        X_err_from_amp_f_sat_dn = zeros(size(mag_X_f));
        X_err_from_amp_s_sat_up = zeros(size(mag_X_s));
        X_err_from_amp_s_sat_dn = zeros(size(mag_X_s));
        X_err_from_amp_u_sat    = zeros(size(mag_X_u));
    end

    constr_viol_err = Inf(size(mag_X_f));
    

    % at sea states where one of the saturated solutions necessarily applies, 
    % set constraint violation error as deviation from that solution,
    % and optimality error as deviation from the optimal ctrl_mult_phase.
    % 
    constr_viol_err(idx_force_sat_applies)    = F_err_from_force_sat(idx_force_sat_applies);
    constr_viol_err(idx_amp_f_sat_applies_up) = X_err_from_amp_f_sat_up(idx_amp_f_sat_applies_up);
    constr_viol_err(idx_amp_f_sat_applies_dn) = X_err_from_amp_f_sat_dn(idx_amp_f_sat_applies_dn);
    constr_viol_err(idx_amp_s_sat_applies_up) = X_err_from_amp_s_sat_up(idx_amp_s_sat_applies_up);
    constr_viol_err(idx_amp_s_sat_applies_dn) = X_err_from_amp_s_sat_dn(idx_amp_s_sat_applies_dn);
    constr_viol_err(idx_amp_u_sat_applies)    = X_err_from_amp_u_sat(idx_amp_u_sat_applies);
    constr_viol_err(idx_power_sat_applies_up) = P_err_from_power_sat_up(idx_power_sat_applies_up);
    constr_viol_err(idx_power_sat_applies_dn) = P_err_from_power_sat_dn(idx_power_sat_applies_dn);
    
    % at sea states where none of the solutions necessarily apply (the true 
    % solution could be any combo of one or two of the saturated solutions but unclear which 1-2),
    % set one error as the avg constraint violation (to penalize exceeding the constraint
    % without penalizing going under the constraint, so does not maximize power)
    % and the other error as the deviation from the max-power point when no
    % constraints are active or zero when a constraint is active (so that
    % when in the feasible region, it walks toward the max-power point
    % until it hits a constraint). To prevent it from starting infeasible
    % and terminating at the nearest constraint without ever going to the 
    % interior power-maximizing part, set optimality error=10 when any constraint is
    % violated. This will make it take larger steps when approaching from 
    % the outside so hopefully it overshoots and re-approaches from the inside.

    if isinf(P_max)
        constr_violation = (max(F_err_from_force_sat,   0) + max(X_err_from_amp_f_sat_up,0) + ...
                            max(X_err_from_amp_f_sat_dn,0) + max(X_err_from_amp_s_sat_up,0) + ...
                            max(X_err_from_amp_s_sat_dn,0) + max(X_err_from_amp_u_sat,   0) + ...
                            max(P_err_from_power_sat_dn))  / 7;
    else
        constr_violation = (max(F_err_from_force_sat,   0) + max(X_err_from_amp_f_sat_up,0) + ...
                            max(X_err_from_amp_f_sat_dn,0) + max(X_err_from_amp_s_sat_up,0) + ...
                            max(X_err_from_amp_s_sat_dn,0) + max(X_err_from_amp_u_sat,   0) + ...
                            max(P_err_from_power_sat_dn,0) + max(P_err_from_power_sat_up,0)   )  / 8;
    end
    constr_viol_err(mult_const_viol_unsat) = constr_violation(mult_const_viol_unsat);
    
    % prevent 0/0=NaN when limit of zero is set
    if F_max==0
        constr_viol_err(mag_U==0) = 0;
        constr_viol_err(mag_U~=0) = mag_U(mag_U~=0) ./ mag_U_unsat(mag_U~=0);
    else
        % add penalty for B_p<0 (negative power)
        tol = 0;
        B_p_violation = max(-B_p_sat+tol,0);
        constr_viol_err = constr_viol_err + B_p_violation;
    end

    if any(~isfinite(constr_viol_err(~isnan(H))),'all')
        warning('error is non finite for non-nan sea state')
    end
    if nargout > 1
        optimality_err  = Inf(size(mag_X_f));

        idx_soln_applies = idx_force_sat_applies    | idx_amp_f_sat_applies_up |...
                           idx_amp_f_sat_applies_dn | idx_amp_s_sat_applies_up |...
                           idx_amp_s_sat_applies_dn | idx_amp_u_sat_applies    | ...
                           idx_power_sat_applies_up | idx_power_sat_applies_dn;
        % eqns below from doi:10.1016/j.ifacol.2024.10.093
        alpha = imag(Z_th) ./ real(Z_th); % eqn 8
        epsilon = zeros(size(mag_X_f)); 
        epsilon(idx_force_sat_applies)  = 1; % effort limit
        epsilon((idx_power_sat_applies_up | idx_power_sat_applies_dn) & ~idx_force_sat_applies) = 1; % power limit phase is undefined, so treat it as an effort limit to minimize force for a given power limit.
        epsilon(idx_soln_applies & ~idx_force_sat_applies & ~(idx_power_sat_applies_up | idx_power_sat_applies_dn)) = -1; % flow limit
        sigma = sqrt( (alpha.^2 .* ctrl_mult_mag.^2 + 1).^2 + alpha.^2 .* (ctrl_mult_mag.^2 + 1).^2);
        acos_argument = -2 * alpha .* ctrl_mult_mag ./ sigma;
        acos_argument = max(-1, min(1, acos_argument));
        atan_denominator = sigma + epsilon .* alpha .* (1+ctrl_mult_mag).^2;
        small_denominator = abs(atan_denominator) < eps;
        atan_denominator(small_denominator) = eps .* sign(atan_denominator(small_denominator));
        atan_denominator(small_denominator & atan_denominator == 0) = eps;
        atan_argument = (alpha.^2 .*  ctrl_mult_mag.^2 + 1) ./ atan_denominator;
        ctrl_mult_phase_desired = 2 * atan(atan_argument) + epsilon .* acos(acos_argument); % eqn 4
        optimality_err(idx_soln_applies) = wrapToPi(ctrl_mult_phase(idx_soln_applies) - ctrl_mult_phase_desired(idx_soln_applies));

        P_unsat_safe = max(P_unsat, eps);
        power_lost = (P_unsat - P_sat) ./ P_unsat_safe;

        force_const_active = ismembertol(mag_U,   4/pi * F_max);
        min_f_const_active = ismembertol(mag_X_f, X_f_lower_limit);
        max_f_const_active = ismembertol(mag_X_f, X_f_upper_limit);
        min_s_const_active = ismembertol(mag_X_s, X_s_lower_limit);
        max_s_const_active = ismembertol(mag_X_s, X_s_upper_limit);
        max_u_const_active = ismembertol(mag_X_u, X_u_upper_limit);
        max_p_const_active = ismembertol(P_sat,   P_max);
        min_p_const_active = ismembertol(P_sat,   P_min);

        any_constr_active = force_const_active | min_f_const_active | max_f_const_active ...
                          | min_s_const_active | max_s_const_active | max_u_const_active ...
                          | max_p_const_active | min_p_const_active;
        any_constr_violated = constr_violation > 0 & ~any_constr_active;
        all_constr_inactive = ~any_constr_active & ~any_constr_violated;
        
        optimality_err(mult_const_viol_unsat & any_constr_active) = 0;
        optimality_err(mult_const_viol_unsat & any_constr_violated) = 10;
        optimality_err(mult_const_viol_unsat & all_constr_inactive) = ...
            power_lost(mult_const_viol_unsat & all_constr_inactive);
        
    end
end

function ang = my_angle(complex_term) % also defined in get_response_drag
    ang = atan2(imag(complex_term),real(complex_term));
end
