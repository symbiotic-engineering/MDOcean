function dLCOEdx = simulation_gradient(X, p, P_avg_elec, capex_design, matl)
%SIMULATION_GRADIENT Compute dLCOE/dx using analytical CAPEX and sparse FD power terms.

    hr_per_yr = 8766;
    alpha_opex = 0.5567;
    alpha_pto = 0.206;

    dcapexdx = zeros(length(X),1);
    coeff_force = (0.0086 + 0.0118 * p.N_WEC^(-alpha_pto)) * p.cost_perN_mult;
    coeff_power = (0.4454 + 0.9099 * p.N_WEC^(-alpha_pto)) * p.cost_perW_mult;
    idx_F_max = 6;
    idx_P_max = 7;
    dcapex_design_dep_dx_Fmax = coeff_force * 1e6;
    dcapex_design_dep_dx_Pmax = coeff_power * 1e5;
    dcapexdx(idx_F_max) = p.N_WEC * dcapex_design_dep_dx_Fmax;
    dcapexdx(idx_P_max) = p.N_WEC * dcapex_design_dep_dx_Pmax;

    structural_idxs = [1:5, 8:12];
    dPdx = zeros(length(X),1);
    non_structural_idxs = setdiff(1:length(X),structural_idxs);
    for idx = non_structural_idxs
        step = max(1e-8, abs(X(idx)) * 1e-4);
        X_perturbed = X(:);
        X_perturbed(idx) = X_perturbed(idx) + step;
        [~, ~, ~, val_perturbed] = simulation([X_perturbed; matl], p);
        dPdx(idx) = (val_perturbed.power_avg - P_avg_elec) / step;
    end

    opex_per_wec = 1.193e6 * p.N_WEC^-alpha_opex;
    opex = opex_per_wec * p.N_WEC;
    capex_non_design_dep = 12.68e6 * p.N_WEC^(-0.741) + 1.24e6;
    capex_total = p.N_WEC * (capex_design + capex_non_design_dep);
    numerator = p.FCR * capex_total + opex;
    denominator = p.N_WEC * P_avg_elec * p.eff_array * hr_per_yr / 1000;

    dnumeratordx = p.FCR * dcapexdx;
    ddenominatordx = p.N_WEC * p.eff_array * hr_per_yr / 1000 .* dPdx;
    dLCOEdx = (dnumeratordx .* denominator - numerator .* ddenominatordx) ./ (denominator.^2);
end
