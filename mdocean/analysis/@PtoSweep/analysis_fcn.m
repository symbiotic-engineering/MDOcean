function intermed_result_struct = analysis_fcn(p, b)
% analysis_fcn  Run the 2-D F_max / P_max sweep for nominal RM3 geometry.
%
% Sweeps maximum powertrain force (F_max) and maximum powertrain power
% (P_max) over a uniform grid while holding all geometric and structural
% design variables fixed at their nominal values.  At each grid point the
% full simulation is evaluated and the resulting average power, design
% cost, LCOE, and nonlinear-constraint vector are stored.
%
% :param p: Parameter struct (from parameters())
% :param b: Design variable bounds struct (from var_bounds())
% :returns: intermed_result_struct  Struct with sweep vectors, result
%           matrices, violated flag, and per-sea-state constraint
%           activity counts (n_force_only, n_power_only, n_both, n_neither).

    n = 20;

    % F_max sweep vector (design-variable units: MN)
    F_max_vec = linspace(0.5, 10, n);
    % P_max sweep vector (design-variable units: 100 kW)
    P_max_vec = linspace(0.5, 10, n);

    [F_MAX, P_MAX] = ndgrid(F_max_vec, P_max_vec);

    % Nominal design vector (geometry and structural DVs fixed)
    X_nom = [b.X_noms; 1];
    idx_F = find(strcmp(b.var_names, 'F_max'));
    idx_P = find(strcmp(b.var_names, 'P_max'));

    power_avg    = zeros(size(F_MAX));
    cost         = zeros(size(F_MAX));
    LCOE_mat     = zeros(size(F_MAX));
    violated     = false(size(F_MAX));
    n_force_only = zeros(size(F_MAX));
    n_power_only = zeros(size(F_MAX));
    n_both       = zeros(size(F_MAX));
    n_neither    = zeros(size(F_MAX));

    for i = 1:numel(F_MAX)
        disp(['Running ' num2str(i) ' of ' num2str(numel(F_MAX))])

        X_i        = X_nom;
        X_i(idx_F) = F_MAX(i);
        X_i(idx_P) = P_MAX(i);

        [~, ~, g_vec, val] = simulation(X_i, p);

        power_avg(i) = val.power_avg;
        cost(i)      = val.capex_design;
        LCOE_mat(i)  = val.LCOE;
        violated(i)  = any(g_vec < 0);

        % --- Per-sea-state constraint activity ---
        % Force constraint active: force saturation is reducing power (P_sat_ratio < 1)
        % Power constraint active: unsaturated electrical power exceeds P_max
        force_sat_tol = 1e-6;                                         % numerical tolerance for force saturation
        P_max_W      = X_i(idx_P) * 1e5;                             % actual P_max in W
        valid_ss     = p.JPD(:) > 0;
        force_active = val.P_sat_ratio(:) < (1 - force_sat_tol) & ~isnan(val.P_sat_ratio(:));
        power_active = val.P_mech(:) * p.eff_pto > P_max_W          & ~isnan(val.P_mech(:));

        n_force_only(i) = sum( force_active & ~power_active & valid_ss);
        n_power_only(i) = sum(~force_active &  power_active & valid_ss);
        n_both(i)       = sum( force_active &  power_active & valid_ss);
        n_neither(i)    = sum(~force_active & ~power_active & valid_ss);
    end

    intermed_result_struct.F_max_vec    = F_max_vec;
    intermed_result_struct.P_max_vec    = P_max_vec;
    intermed_result_struct.F_MAX        = F_MAX;
    intermed_result_struct.P_MAX        = P_MAX;
    intermed_result_struct.power_avg    = power_avg;
    intermed_result_struct.cost         = cost;
    intermed_result_struct.LCOE         = LCOE_mat;
    intermed_result_struct.violated     = violated;
    intermed_result_struct.n_force_only = n_force_only;
    intermed_result_struct.n_power_only = n_power_only;
    intermed_result_struct.n_both       = n_both;
    intermed_result_struct.n_neither    = n_neither;
end
