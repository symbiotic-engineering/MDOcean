function intermed_result_struct = analysis_fcn(p, b)
% analysis_fcn  Run the geometry sweep simulation loop.
%
% Sweeps six geometric ratios (a_1/a_2, a_2/h, d_1/h, d_2/d_1, a_3/a_1)
% and water depth h, computing the hydro ratio CW/CW_max and the wave
% parameter m_0*h at every grid point.
%
% :param p: Parameter struct (from parameters())
% :param b: Design variable bounds struct (from var_bounds())
% :returns: intermed_result_struct  Struct with hydro_ratio_result,
%           m0h_stored, val (simulation outputs), all grid arrays, the
%           sweep parameter vectors, and p.

    warning_ids = { ...
        'MDOcean:GetResponseDrag:OpenLoopUnstable', ...
        'MDOcean:RunMEEM:NegativeAddedMass', ...
        'MDOcean:RunMEEM:NegativeDamping', ...
        'MDOcean:RunMEEM:NegativeDampingGlobal', ...
        'MDOcean:RunMEEM:InvalidGeometry', ...
        'MDOcean:RunMEEM:NonFiniteA', ...
        'MDOcean:RunMEEM:NonFiniteB', ...
        'MDOcean:RunMEEM:NonFiniteC' ...
    };

    warning_ids_with_state = [{'MATLAB:nearlySingularMatrix', 'MATLAB:singularMatrix', 'backtrace'}, warning_ids];
    warning_states = cell(size(warning_ids_with_state));
    for k = 1:numel(warning_ids_with_state)
        warning_states{k} = warning('query', warning_ids_with_state{k});
    end
    cleanup_warnings = onCleanup(@() restore_warning_states(warning_ids_with_state, warning_states));

    warning('off','MATLAB:nearlySingularMatrix')
    warning('off','MATLAB:singularMatrix')
    warning('off','backtrace') % disable longer warning messages
    for k = 1:numel(warning_ids)
        warning('error', warning_ids{k})
    end

    n = 3; % number of points per dimension

    zero_to_one   = linspace(0.01, 0.99, n);
    zero_to_three = linspace(0.01, 3,    n);

    h_vec = [10 200]; % gives a range of m0h of 0.36-2.06 and 2.40-39.79 for default wave periods
    a1_a2 = zero_to_one;
    a2_h  = zero_to_three;
    d1_h  = zero_to_one;
    d2_d1 = zero_to_one;
    a3_a1 = 1 + zero_to_three;

    [H, A1_A2, A2_H, D1_H, D2_D1, A3_A1] = ndgrid(h_vec, a1_a2, a2_h, d1_h, d2_d1, a3_a1);

    A1_H = A1_A2 .* A2_H;
    D2_H = D2_D1 .* D1_H;
    A3_H = A3_A1 .* A1_H;

    A1 = A1_H .* H;
    A2 = A2_H .* H;
    A3 = A3_H .* H;
    D1 = D1_H .* H;
    D2 = D2_H .* H;

    X = [b.X_noms; 1];

    nT = length(p.T);
    n_geoms = numel(A1);
    m0h_stored_linear         = zeros(nT, n_geoms);
    hydro_ratio_result_linear = nan(nT, n_geoms);
    warning_hits              = false(numel(warning_ids), n_geoms);
    val = repmat(struct('vol_f', NaN, 'vol_s', NaN), [1, n_geoms]);

    parfor i = 1:n_geoms

        D_s   = 2 * A1(i);
        D_f   = 2 * A2(i);
        D_d   = 2 * A3(i);
        T_s   = D1(i);
        T_f_2 = D2(i);

        % Construct per-iteration parameter struct and design vector
        p_i = p;
        p_i.h              = H(i);
        p_i.T_s_over_D_s   = T_s   / D_s;
        p_i.D_d_over_D_s   = D_d   / D_s;
        X_i = X;
        X_i(strcmp(b.var_names, 'D_s'))   = D_s;
        X_i(strcmp(b.var_names, 'D_f'))   = D_f;
        X_i(strcmp(b.var_names, 'h_s'))   = T_s + 5;
        X_i(strcmp(b.var_names, 'T_f_2')) = T_f_2;

        m0h_stored_linear(:, i) = dispersion(2*pi ./ p.T, p_i.h, p.g) * p_i.h;

        [hydro_ratio_max, out, warning_hit] = run_check_max_CW_with_warning_capture(p_i, X_i, warning_ids);
        val(i) = ensure_vol_fields(out);
        warning_hits(:, i) = warning_hit(:);
        hydro_ratio_result_linear(:, i) = hydro_ratio_max;
    end

    m0h_stored = reshape(m0h_stored_linear, [nT, size(A1)]);
    hydro_ratio_result = reshape(hydro_ratio_result_linear, [nT, size(A1)]);

    restore_warning_states(warning_ids_with_state, warning_states)
    clear cleanup_warnings

    warning_counts = sum(warning_hits, 2);
    for k = 1:numel(warning_ids)
        if warning_counts(k) > 0
            warning('MDOcean:SweepGeoms:WarningSummary', ...
                'SweepGeoms caught %s for %.1f%% of geometries (%d/%d).', ...
                warning_ids{k}, 100 * warning_counts(k) / n_geoms, warning_counts(k), n_geoms)
        end
    end

    % --- Save intermediate results ---
    % Main simulation outputs
    intermed_result_struct.hydro_ratio_result = hydro_ratio_result;
    intermed_result_struct.m0h_stored         = m0h_stored;
    intermed_result_struct.val                = val;

    % Dimensional grid arrays (needed for CWR computation)
    intermed_result_struct.H     = H;
    intermed_result_struct.A1_A2 = A1_A2;
    intermed_result_struct.A2_H  = A2_H;
    intermed_result_struct.D1_H  = D1_H;
    intermed_result_struct.D2_D1 = D2_D1;
    intermed_result_struct.A3_A1 = A3_A1;
    intermed_result_struct.A1    = A1;
    intermed_result_struct.A2    = A2;
    intermed_result_struct.D1    = D1;
    intermed_result_struct.D2    = D2;

    % Sweep parameter vectors (needed for plot annotations)
    intermed_result_struct.h_vec  = h_vec;
    intermed_result_struct.a1_a2  = a1_a2;
    intermed_result_struct.a2_h   = a2_h;
    intermed_result_struct.d1_h   = d1_h;
    intermed_result_struct.d2_d1  = d2_d1;
    intermed_result_struct.a3_a1  = a3_a1;

    % Parameter struct (p.T and p.g needed for CWR in post-processing)
    intermed_result_struct.p = p;
end

function [hydro_ratio_max, out, warning_hit] = run_check_max_CW_with_warning_capture(p_i, X_i, warning_ids)
    warning_hit = false(1, numel(warning_ids));
    disabled_warning_ids = {};

    while true
        try
            [hydro_ratio, ~, ~, ~, ~, ~, ~, out] = check_max_CW('', p_i, X_i, false);
            hydro_ratio_max = max(hydro_ratio);
            break
        catch ME
            warning_idx = find(strcmp(ME.identifier, warning_ids), 1);
            if isempty(warning_idx)
                rethrow(ME)
            end
            warning_hit(warning_idx) = true;

            if any(strcmp(ME.identifier, disabled_warning_ids))
                rethrow(ME)
            end
            warning('off', ME.identifier)
            disabled_warning_ids{end + 1} = ME.identifier; %#ok<AGROW>
        end
    end

    for k = 1:numel(disabled_warning_ids)
        warning('error', disabled_warning_ids{k})
    end
end

function out = ensure_vol_fields(out)
    if ~isfield(out, 'vol_f')
        out.vol_f = NaN;
    end
    if ~isfield(out, 'vol_s')
        out.vol_s = NaN;
    end
end

function restore_warning_states(warning_ids, warning_states)
    for k = 1:numel(warning_ids)
        warning(warning_states{k}.state, warning_ids{k})
    end
end
