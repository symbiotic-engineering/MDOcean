function [CEM_CO2, CEM_wec_capacity, CEM_grid_cost] = grid_CEM(B_p, X_u, phase_X_u, ...
                                                gamma_phase_f, gamma_f_over_rho_g, capacity_cost, location)


    [zeta, omega_n] = fit_second_order_sys(X_u, phase_X_u, gamma_f_over_rho_g, gamma_phase_f);

    [CEM_CO2, CEM_wec_capacity, CEM_grid_cost] = CEM_lookup_table(zeta, omega_n, capacity_cost, B_p, location);

end


function [zeta, omega_n] = fit_second_order_sys(X_u, phase_X_u, gamma_f_over_rho_g, gamma_phase_f)
    % fixme put real fit here

    %combined()

    %[zeta, omega_n] = combined(X_u, phase_X_u, gamma_phase_f) don't use
    %combined.m here

    
    
    [zeta, omega_n] = fit_from_vars(X_u, phase_X_u, gamma_f_over_rho_g, gamma_phase_f);

    %zeta
    %omega_n

    %zeta = 0.05;
    %omega_n = 0.4;

end

function row = findNearestRow(zeta0, omega_n0, capCost0, T)


    % Extract numeric arrays
    z = T.zeta;
    w = T.omega_n;
    c = T.capacity_cost;

    % Check ranges
    if zeta0 < min(z) || zeta0 > max(z)
       disp(zeta0)
        error('zeta inputs are out of range.');
        
    end

    if   omega_n0 < min(w) || omega_n0 > max(w)
        error('omega are out of range.');
    end

    if capCost0 < min(c) || capCost0 > max(c)
        error('capCost out of range')

    end

    % find distance (placeholder interp code, fixme)
    difference = (z - zeta0).^2 + (w - omega_n0).^2 + (c - capCost0).^2;

    % Find index of minimum distance
    [~, idx] = min(difference);

    % Return that row
    row = T(idx, :);
end



function [CEM_CO2, CEM_wec_capacity, CEM_grid_cost] = CEM_lookup_table(zeta, omega_n, capacity_cost, B_p, location)

    dummy = readtable('dummy3.csv') % this is just placeholder for now, replace with actual data
    %{
    CEM_co2 = interpn(zeta_data, omega_n_data, cap_cost_data, power_lim_data, CO2_data,...
        zeta,omega_n,capacity_cost,power_limit);
    CEM_capacity = interpn(zeta_data, omega_n_data, cap_cost_data, power_lim_data, capacity_data,...
        zeta,omega_n,capacity_cost,power_limit);
    %}

    %CEM_co2 = interpn(dummy.zeta, dummy.omega_n, cap_cost_data, power_lim_data, CO2_data,...
    %    zeta,omega_n,capacity_cost,power_limit);
    %CEM_capacity = interpn(zeta_data, omega_n_data, cap_cost_data, power_lim_data, capacity_data,...
    %    zeta,omega_n,capacity_cost,power_limit);

    
    %if zeta == 0.05 && omega_n == 0.4 && strcmp(location,'ISONE')
    
    CEM_data = findNearestRow(zeta, omega_n, capacity_cost, dummy)

    disp(CEM_data)

    isValidLookupLocation = true;% placeholder
    
    if (isValidLookupLocation)
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