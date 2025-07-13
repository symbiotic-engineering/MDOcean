function [CEM_CO2, CEM_wec_capacity, CEM_grid_cost] = grid_CEM(B_p, X_u, phase_X_u, ...
                                                gamma_phase_f, capacity_cost, location)


    [zeta, omega_n] = fit_second_order_sys(X_u, phase_X_u, gamma_phase_f);

    [CEM_CO2, CEM_wec_capacity, CEM_grid_cost] = CEM_lookup_table(zeta, omega_n, capacity_cost, B_p, location);

end


function [zeta, omega_n] = fit_second_order_sys(X_u, phase_X_u, gamma_phase_f)
    % fixme put real fit here
    zeta = 0.05;
    omega_n = 0.4;

end

function [CEM_CO2, CEM_wec_capacity, CEM_grid_cost] = CEM_lookup_table(zeta, omega_n, capacity_cost, B_p, location)
    if zeta == 0.05 && omega_n == 0.4 && strcmp(location,'ISONE')
        % fixme put real lookup table here
        CEM_CO2 = 1 + capacity_cost;
        CEM_wec_capacity = 1 - capacity_cost;
        CEM_grid_cost = 1 - capacity_cost;
    else
        error('CEM results not available for this zeta, w_n, and location')
    end

end