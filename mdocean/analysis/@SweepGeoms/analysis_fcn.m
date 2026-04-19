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
        'MDOcean:RunMEEM:NonFiniteC', ...
        'MDOcean:DragIntegral:OutsideLUT', ...
        'MDOcean:GetResponseDrag:StabilizingFailed', ...
        'MDOcean:Simulation:InfNanImaginary', ...
        'MATLAB:integral:NonFiniteValue' ...
    };

    verbosity_warning_ids = {'MATLAB:nearlySingularMatrix', 'MATLAB:singularMatrix', 'backtrace'};

    % Save all current warning states and restore on function exit.
    prev_warning_states = warning;
    cleanup_warnings = onCleanup(@() warning(prev_warning_states)); %#ok<NASGU>

    % Suppress noisy verbosity warnings in the outer scope (nominal simulation
    % and post-loop code).  Per-geometry warnings inside the parfor are handled
    % silently by run_and_catch_warnings via evalc and do not need to be set
    % here.
    for k = 1:numel(verbosity_warning_ids)
        warning('off', verbosity_warning_ids{k});
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

    % Run one nominal simulation to discover the full output struct layout
    % for pre-allocating the val array before the parfor loop.
    [~, ~, ~, val_nom] = simulation(X, p);
    fn = fieldnames(val_nom);
    nan_val = val_nom;
    for k = 1:numel(fn)
        nan_val.(fn{k}) = NaN(size(val_nom.(fn{k})));
    end

    nT = length(p.T);
    n_geoms = numel(A1);
    m0h_stored_linear         = nan(nT, n_geoms);
    hydro_ratio_result_linear = nan(nT, n_geoms);
    warning_hits              = false(numel(warning_ids), n_geoms);
    unknown_error_hits        = false(1, n_geoms);
    val = repmat(nan_val, [1, n_geoms]);
    drag_lut_rp_range    = nan(2, n_geoms); % row 1 = min, row 2 = max rp per geometry
    drag_lut_kappa_range = nan(2, n_geoms); % row 1 = min, row 2 = max kappa per geometry

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

        try
            [hydro_ratio_max, out, warning_hit, rp_range_i, kappa_range_i] = run_check_max_CW_with_warning_capture(p_i, X_i, warning_ids, nan_val);
        catch
            hydro_ratio_max  = NaN;
            out              = nan_val;
            warning_hit      = false(1, numel(warning_ids));
            rp_range_i       = nan(2, 1);
            kappa_range_i    = nan(2, 1);
            unknown_error_hits(i) = true;
        end
        val(i) = out;
        warning_hits(:, i) = warning_hit(:);
        hydro_ratio_result_linear(:, i) = hydro_ratio_max;
        drag_lut_rp_range(:, i)    = rp_range_i;
        drag_lut_kappa_range(:, i) = kappa_range_i;
    end

    m0h_stored = reshape(m0h_stored_linear, [nT, size(A1)]);
    hydro_ratio_result = reshape(hydro_ratio_result_linear, [nT, size(A1)]);

    warning_counts = sum(warning_hits, 2);
    drag_lut_id = 'MDOcean:DragIntegral:OutsideLUT';
    for k = 1:numel(warning_ids)
        if warning_counts(k) > 0
            if strcmp(warning_ids{k}, drag_lut_id)
                valid_cols  = ~isnan(drag_lut_rp_range(1, :));
                rp_min_all    = min(drag_lut_rp_range(1, valid_cols));
                rp_max_all    = max(drag_lut_rp_range(2, valid_cols));
                kappa_min_all = min(drag_lut_kappa_range(1, valid_cols));
                kappa_max_all = max(drag_lut_kappa_range(2, valid_cols));
                warning('MDOcean:SweepGeoms:WarningSummary', ...
                    ['SweepGeoms caught %s for %.1f%% of geometries (%d/%d). ' ...
                     'rp range across all geometries: [%.6g, %.6g], ' ...
                     'kappa range: [%.6g, %.6g].'], ...
                    warning_ids{k}, 100 * warning_counts(k) / n_geoms, warning_counts(k), n_geoms, ...
                    rp_min_all, rp_max_all, kappa_min_all, kappa_max_all)
            else
                warning('MDOcean:SweepGeoms:WarningSummary', ...
                    'SweepGeoms caught %s for %.1f%% of geometries (%d/%d).', ...
                    warning_ids{k}, 100 * warning_counts(k) / n_geoms, warning_counts(k), n_geoms)
            end
        end
    end
    n_unknown_errors = sum(unknown_error_hits);
    if n_unknown_errors > 0
        warning('MDOcean:SweepGeoms:UnknownError', ...
            'SweepGeoms silently suppressed unknown errors for %.1f%% of geometries (%d/%d).', ...
            100 * n_unknown_errors / n_geoms, n_unknown_errors, n_geoms);
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

function [hydro_ratio_max, out, warning_hit, rp_range, kappa_range] = run_check_max_CW_with_warning_capture(p_i, X_i, warning_ids, nan_val)
    hydro_ratio_max = NaN;
    out = nan_val; % default for early-exit paths (genuine non-warning errors)
    rp_range    = nan(2, 1); % [min; max] rp  seen in MDOcean:DragIntegral:OutsideLUT
    kappa_range = nan(2, 1); % [min; max] kappa seen

    % run_and_catch_warnings runs check_max_CW exactly once inside evalc,
    % capturing all Command-Window output (including warning messages) so
    % nothing reaches the terminal.  Warning states are managed internally.
    [check_outputs, warning_hit, captured_text] = run_and_catch_warnings( ...
        @() check_max_CW('', p_i, X_i, false), warning_ids, 8);
    hydro_ratio     = check_outputs{1};
    out             = check_outputs{8};
    hydro_ratio_max = max(hydro_ratio);

    % Parse rp/kappa ranges from the captured warning text for the drag-LUT
    % warning.  The drag-LUT warning may fire multiple times within a single
    % check_max_CW call (for different rp/kappa values), so collect all
    % occurrences and aggregate min/max.
    drag_lut_id = 'MDOcean:DragIntegral:OutsideLUT';
    if contains(captured_text, drag_lut_id)
        all_toks = regexp(captured_text, ...
            'rp_min=([\d.e+\-]+),\s*rp_max=([\d.e+\-]+),\s*kappa_min=([\d.e+\-]+),\s*kappa_max=([\d.e+\-]+)', ...
            'tokens');
        for m = 1:numel(all_toks)
            vals = cellfun(@str2double, all_toks{m});
            rp_range(1)    = min([rp_range(1);    vals(1)], [], 'omitnan');
            rp_range(2)    = max([rp_range(2);    vals(2)], [], 'omitnan');
            kappa_range(1) = min([kappa_range(1); vals(3)], [], 'omitnan');
            kappa_range(2) = max([kappa_range(2); vals(4)], [], 'omitnan');
        end
    end
end
