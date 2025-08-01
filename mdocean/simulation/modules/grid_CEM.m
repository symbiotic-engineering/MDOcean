function [CEM_CO2, CEM_wec_capacity, CEM_grid_cost] = grid_CEM(B_p, X_u, phase_X_u, ...
                                                gamma_phase_f, capacity_cost, location)


    [zeta, omega_n] = fit_second_order_sys(X_u, phase_X_u, gamma_phase_f);

    [CEM_CO2, CEM_wec_capacity, CEM_grid_cost] = CEM_lookup_table(zeta, omega_n, capacity_cost, B_p, location);

end


function [zeta, omega_n] = fit_second_order_sys(X_u, phase_X_u, gamma_phase_f)
    % fixme put real fit here

    %combined()

    [zeta, omega_n] = combined(X_u, phase_X_u, gamma_phase_f)

    %zeta = 0.05;
    %omega_n = 0.4;

end

function [CEM_CO2, CEM_wec_capacity, CEM_grid_cost] = CEM_lookup_table(zeta, omega_n, capacity_cost, B_p, location)
    %CEM_co2 = interpn(zeta_data, omega_n_data, cap_cost_data, power_lim_data, CO2_data,...
    %    zeta,omega_n,capacity_cost,power_limit);
    %CEM_capacity = interpn(zeta_data, omega_n_data, cap_cost_data, power_lim_data, capacity_data,...
    %    zeta,omega_n,capacity_cost,power_limit);
    if zeta == 0.05 && omega_n == 0.4 && strcmp(location,'ISONE')
        % fixme put real lookup table here
        co2_slope = 2/3;
        cost_slope = 1/3;
        cap_slope = 2.5/3;

        cutin_capacity_cost = 750e3;
        cheapest_cost_with_data = 400e3;

        no_wec_CO2 = 12.003064e6; % tonnes (typical value 10e6=10 MT)
        no_wec_grid_cost = 2.889863123e9;
        
        if capacity_cost > cutin_capacity_cost
            % no WECs
            CEM_CO2 = no_wec_CO2;
            CEM_wec_capacity = 0;
            CEM_grid_cost = no_wec_grid_cost;
        elseif capacity_cost > cheapest_cost_with_data
            % some wecs, and in bounds of model
            capacity_cost_pct_incr = (capacity_cost - 725e3) / 725e3;
    
            CEM_CO2 =  7.551749e6 * (1 + capacity_cost_pct_incr * co2_slope);
            CEM_wec_capacity = 10.201e3 * (1 - capacity_cost_pct_incr * cap_slope);
            CEM_grid_cost = 2.468040544e9 * (1 - capacity_cost_pct_incr * cost_slope);
        else 
            % not in bounds of model
            error('WEC is too cheap, no CEM data here.')
        end
    else
        error('CEM results not available for this zeta, w_n, and location')
    end

end