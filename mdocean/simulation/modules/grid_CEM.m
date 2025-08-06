function [CEM_CO2, CEM_wec_capacity, CEM_grid_cost] = grid_CEM(B_p, X_u, phase_X_u, ...
                                                gamma_phase_f, gamma_f_over_rho_g, capacity_cost, location)
% capacity cost is in $/kW = $k/MW


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

function row = findNearestRow_dummy(zeta0, omega_n0, capCost0, T)

    % sort vectors for queried values
    T = sortrows(T, {'zeta','omega_n','capacity_cost'});
    zVals = unique(T.zeta);
    wVals = unique(T.omega_n);
    cVals = unique(T.capacity_cost);

    % Check input ranges
    if zeta0 < min(zVals) || zeta0 > max(zVals)
       disp(zeta0)
        error('zeta inputs are out of range.');
        
    end

    if   omega_n0 < min(wVals) || omega_n0 > max(wVals)
        error('omega are out of range.');
    end

    if capCost0 < min(cVals) || capCost0 > max(cVals)
        error('capCost out of range')

    end

    nZ = numel(zVals);
    nW = numel(wVals);
    nC = numel(cVals);

    co2Grid      = reshape(T.co2capacity, [nZ, nW, nC]);
    gridcostGrid = reshape(T.gridcost, [nZ, nW, nC]);

    % interpolation
    co2_out      = interpn(zVals, wVals, cVals, co2Grid,      zeta0, omega_n0, capCost0, 'linear');
    gridcost_out = interpn(zVals, wVals, cVals, gridcostGrid, zeta0, omega_n0, capCost0, 'linear');
    % wec_capacity_out = interpn(zVals, wVals, cVals, wecGrid, zeta0, omega_n0, capCost0, 'linear');
    row = [zeta0, omega_n0, capCost0, co2_out, gridcost_out];
end
function row = findNearestRow_interp(zeta0, omega_n0, wecCost0, T)


    %check for holes in data
    %{
    mask = (T.power_lim == powerLim0) & (T.diameter == diameter0);
    Tsub = T(mask, :);
    if isempty(Tsub)
        error('No data for power_lim=%.1f & diameter=%.1f.', powerLim0, diameter0);
    end
    %}

    % Extract grid vectors (replace with Tsub when hole code is working)
    zVals = unique(T.zeta);
    wVals = unique(T.omega_n);
    cVals = unique(T.wec_cost);

    % check for bounds
    if zeta0 < min(zVals) || zeta0 > max(zVals)
       disp(zeta0)
        error('zeta inputs are out of range.');
        
    end

    if   omega_n0 < min(wVals) || omega_n0 > max(wVals)
        error('omega are out of range.');
    end

    if wecCost0 < min(cVals) || wecCost0 > max(cVals)
        wecCost0
        error('weccost out of range')

    end

    %sort arrays
    T = sortrows(T, {'zeta','omega_n','wec_cost'});
    nZ = numel(zVals);
    nW = numel(wVals);
    nC = numel(cVals);

    % identify which columns to interpolate:
    allNames = T.Properties.VariableNames;
    gridAxes  = {'zeta','omega_n','wec_cost'};
    fixedCols = {'power_lim','diameter'};
    % Some columns are categorical or identifiers—exclude those
    exclude = [gridAxes fixedCols {'electrification','carbon_constraint','year'}];
    interpVars = setdiff(allNames, exclude, 'stable');

    % create arrays for results
    results = struct();
    % include the inputs in the output row
    %results.power_lim = powerLim0;
    %results.diameter  = diameter0;
    results.zeta      = zeta0;
    results.omega_n   = omega_n0;
    results.wec_cost  = wecCost0;

    for i = 1:numel(interpVars)
        var = interpVars{i};
        % Reshape into [nZ × nW × nC]
        grid3D = reshape(T.(var), [nZ, nW, nC]);
        % Trilinear interpolation
        results.(var) = interpn( ...
            zVals, wVals, cVals, ...
            grid3D, ...
            zeta0, omega_n0, wecCost0, ...
            'linear' ...
        );
    end

    % return row of lookup table results
    row = struct2table(results);


    % write own interpolator?
end

function [CEM_CO2, CEM_wec_capacity, CEM_grid_cost] = CEM_lookup_table(zeta, omega_n, capacity_cost, B_p, location)

    dummy = readtable('dummy3.csv'); % this is just placeholder for now, replace with actual data
    data_V1 = readtable('scenario_outputs.csv');
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
    
    %CEM_data = findNearestRow_interp(zeta, omega_n, capacity_cost, data_V1)

    %disp(CEM_data)

    isValidLookupLocation = true;% placeholder
    
    if (isValidLookupLocation)
        % fixme put real lookup table here


        co2_slope = 2/3;
        cost_slope = 1/3;
        cap_slope = 2.5/3;


        % placeholders
        cutin_capacity_cost = 750e3;
        cheapest_cost_with_data = 400e3;

        no_wec_CO2 = 12.003064e6; % tonnes (typical value 10e6=10 MT)
        no_wec_grid_cost = 2.889863123e9;

        % data



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
            cheapest_cost_with_data
            capacity_cost
            error('WEC is too cheap, no CEM data here.')
        end
    else
        error('CEM results not available for this zeta, w_n, and location')
    end

end