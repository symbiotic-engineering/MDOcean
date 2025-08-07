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
    % combined.m here

    
    
    [zeta, omega_n] = fit_from_vars(X_u, phase_X_u, gamma_f_over_rho_g, gamma_phase_f);

    zeta
    omega_n
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

function row = findNearestRow_interp(zeta0, omega_n0, wecCost0, powerLim0, T)

    % use specificed power limit (fixme, add to interpolation later)
    Tsub = T(T.power_lim == powerLim0, :);
    if isempty(Tsub)
        error('No data for power_lim = %g.', powerLim0);
    end

    % grid axes for interpolated values
    zVals = unique(Tsub.zeta);
    wVals = unique(Tsub.omega_n);
    cVals = unique(Tsub.wec_cost);

    % check variable bounds
    if any([zeta0<min(zVals), zeta0>max(zVals), ...
            omega_n0<min(wVals), omega_n0>max(wVals), ...
            wecCost0<min(cVals), wecCost0>max(cVals)])
        warning('Query outside data bounds.');
    end

    % sort values
    Tsub = sortrows(Tsub, {'zeta','omega_n','wec_cost'});
    nZ = numel(zVals);
    nW = numel(wVals);
    nC = numel(cVals);

    % interpolated values
    %varsAll   = Tsub.Properties.VariableNames;
    %gridAxes  = {'power_lim','zeta','omega_n','wec_cost'};
    %interpVars = setdiff(varsAll, gridAxes, 'stable');
    %interpVars = 


    % 5a) Grab only the numeric variables
    numTable   = Tsub(:, vartype('numeric'));  
    numNames   = numTable.Properties.VariableNames;

    % 5b) Exclude the three grid axes
    gridAxes   = {'zeta','omega_n','wec_cost'};      
    interpVars = setdiff(numNames, gridAxes, 'stable');

    % get interpolated values
    S = struct( ...
      'power_lim', powerLim0, ...
      'zeta',      zeta0, ...
      'omega_n',   omega_n0, ...
      'wec_cost',  wecCost0);

    for i = 1:numel(interpVars)
        v = interpVars{i};
        % reshape into [nZ × nW × nC]
        G = reshape(Tsub.(v), [nZ, nW, nC]);
        S.(v) = interpn( ...
            zVals, wVals, cVals, ...
            G, ...
            zeta0, omega_n0, wecCost0, ...
            'spline' ...
        );
    end

    % return results
    row = struct2table(S);
end


function [CEM_CO2, CEM_wec_capacity, CEM_grid_cost] = CEM_lookup_table(zeta, omega_n, capacity_cost, B_p, location)

    data_V1 = readtable('scenario_outputs.csv');
    data_V1.wec_cost( data_V1.wec_cost == 5000 ) = 15000; % fixme: this is just placeholder for testing

    %{
    CEM_co2 = interpn(zeta_data, omega_n_data, cap_cost_data, power_lim_data, CO2_data,...
        zeta,omega_n,capacity_cost,power_limit);
    CEM_capacity = interpn(zeta_data, omega_n_data, cap_cost_data, power_lim_data, capacity_data,...
        zeta,omega_n,capacity_cost,power_limit);
    %}

    
    %if zeta == 0.05 && omega_n == 0.4 && strcmp(location,'ISONE')
    
    powerLim_dummy = 700.0 % fixme: add to simulation later on

    zeta
    omega_n
    capacity_cost

    CEM_data = findNearestRow_interp(zeta, omega_n, capacity_cost, powerLim_dummy, data_V1)

    isValidLookupLocation = true;% placeholder
    
    if (isValidLookupLocation)
        % fixme put real lookup table here


        %co2_slope = 2/3;
        %cost_slope = 1/3;
        %cap_slope = 2.5/3;


        % placeholders
        cutin_capacity_cost = 750;
        cheapest_cost_with_data = 400;

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
            %capacity_cost_pct_incr = (capacity_cost - 725e3) / 725e3;
    
            %CEM_CO2 =  7.551749e6 * (1 + capacity_cost_pct_incr * co2_slope);
            %CEM_wec_capacity = 10.201e3 * (1 - capacity_cost_pct_incr * cap_slope);
            %CEM_grid_cost = 2.468040544e9 * (1 - capacity_cost_pct_incr * cost_slope);

            CEM_CO2 = CEM_data.carbon;
            CEM_wec_capacity = CEM_data.wave_capacity;
            CEM_grid_cost = CEM_data.system_cost;

            
        else 
            % not in bounds of model

            %cheapest_cost_with_data
            capacity_cost
            error('WEC is too cheap, no CEM data here.')
        end
    else
        error('CEM results not available for this zeta, w_n, and location')
    end

end