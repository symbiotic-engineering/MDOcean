function [CEM_CO2, CEM_wec_capacity, CEM_grid_cost] = grid_CEM(B_p, X_u, phase_X_u, ...
                                                gamma_phase_f, gamma_f_over_rho_g, capacity_cost, power_lim_frac, location, params)
% capacity cost is in $/kW = $k/MW


    [zeta, omega_n] = fit_second_order_sys(X_u, phase_X_u, gamma_f_over_rho_g, gamma_phase_f, params);

    [CEM_CO2, CEM_wec_capacity, CEM_grid_cost] = CEM_lookup_table(zeta, omega_n, capacity_cost, power_lim_frac, B_p, location, params);

end


function [zeta, omega_n] = fit_second_order_sys(X_u, phase_X_u, gamma_f_over_rho_g, gamma_phase_f, params)
    % fixme put real fit here

    %combined()

    %[zeta, omega_n] = combined(X_u, phase_X_u, gamma_phase_f) don't use
    % combined.m here

    
    
    [zeta, omega_n] = fit_from_vars(X_u, phase_X_u, gamma_f_over_rho_g, gamma_phase_f, params);


    %zeta = 0.05;
    %omega_n = 0.4;


end



function row = findNearestRow_interp(zeta0, omega_n0, wecCost0, powerLim0, T)

    % grid axes for interpolated values
    zVals = unique(T.zeta);
    wVals = unique(T.omega_n);
    cVals = unique(T.wec_cost);
    pVals = unique(T.power_lim);

    % check variable bounds
    if any([zeta0<min(zVals), zeta0>max(zVals), ...
            omega_n0<min(wVals), omega_n0>max(wVals), ...
            wecCost0<min(cVals), wecCost0>max(cVals),...
            powerLim0<min(pVals), powerLim0>max(pVals)])
        fprintf(['Query outside data bounds, with values ' ...
            'zeta=%.3f, omega_n=%.3f, cost=%.3f, power lim=%.3f\n'], ...
            zeta0, omega_n0, wecCost0, powerLim0);
    end

    % sort values
    gridAxes   = {'zeta','omega_n','wec_cost','power_lim'};
    T = sortrows(T, gridAxes);
    nZ = numel(zVals);
    nW = numel(wVals);
    nC = numel(cVals);
    nP = numel(pVals);

    % interpolated values
    %varsAll   = Tsub.Properties.VariableNames;
    %gridAxes  = {'power_lim','zeta','omega_n','wec_cost'};
    %interpVars = setdiff(varsAll, gridAxes, 'stable');
    %interpVars = 


    % 5a) Grab only the numeric variables
    numTable   = T(:, vartype('numeric'));  
    numNames   = numTable.Properties.VariableNames;

    % 5b) Exclude the three grid axes
    interpVars = setdiff(numNames, gridAxes, 'stable');

    % get interpolated values
    S = struct( ...
      'power_lim', powerLim0, ...
      'zeta',      zeta0, ...
      'omega_n',   omega_n0, ...
      'wec_cost',  wecCost0);

    idx_not_nan = ~isnan(T.wave_capacity);

    for i = 1:numel(interpVars)
        v = interpVars{i};
        % reshape into [nZ × nW × nC]
        % this method cannot extrapolate, any dimension number ok, requires
        % full grid.
        G = reshape(T.(v), [nZ, nW, nC, nP]);
%         S.(v) = interpn( ...
%             zVals, wVals, cVals, ...
%             G, ...
%             zeta0, omega_n0, wecCost0, ...
%             'spline' ...
%         );

% this method is good for when there are nans, so not a full grid. can extrapolate.
% but can only handle 3 input dimensions (so not power limit)
%         interpFcn = scatteredInterpolant(T.zeta(idx_not_nan), ...
%                                         T.omega_n(idx_not_nan), ...
%                                         T.wec_cost(idx_not_nan), ...
%                                         T.(v)(idx_not_nan));
%         S.(v) = interpFcn(zeta0, omega_n0, wecCost0);

% this method requires a full grid and can extrapolate. any dimensions # ok.
% janky assumption: nans are zero wec capacity
        G_zero_capacity = G(T.wave_capacity==0);
        G_filled_nan = G;
        G_filled_nan(isnan(G)) = G_zero_capacity(1);
        interpFcn = griddedInterpolant({zVals, ...
                                       wVals, ...
                                       cVals, ...
                                       pVals}, ...
                                       G_filled_nan);
        S.(v) = interpFcn(zeta0, omega_n0, wecCost0, powerLim0);

    end

    % return results
    row = struct2table(S);

    unviable = row.wave_capacity<=0;
    bad_extrap = any(row.Variables<0); % negative numbers meaningless
    %{
    if unviable || bad_extrap
        % if not viable, use margin to viability instead (how much cost
        % needs to decrease in order to be viable)

        idx_viable = T.wave_capacity > 0;

        % scattered interpolant only takes 3D, so limit to nearest

        capacity_viable = T.wave_capacity(idx_viable);
        zeta_viable     = T.zeta(idx_viable);
        omega_viable    = T.omega_n(idx_viable);
        cost_viable     = T.wec_cost(idx_viable);

        try
            costViabilityThresholdFcn = scatteredInterpolant(zeta_viable, omega_viable, capacity_viable, cost_viable);
        catch ME
            if (strcmp(ME.identifier,'MATLAB:griddedInterpolant:DegenerateGridErrId'))
                msg = 'There is not enough CEM data with nonzero capacity.';
                causeException = MException('MATLAB:MDOcean:CEM_nonzero',msg);
                ME = addCause(ME,causeException);
            end
            rethrow(ME)
        end
        wecCostThresholdViable = costViabilityThresholdFcn(zeta0, omega_n0, 0);
        margin_to_viability = wecCost0 - wecCostThresholdViable;
        assert(margin_to_viability>0)
    end
    %}
end


function [CEM_CO2, CEM_wec_capacity, CEM_grid_cost] = CEM_lookup_table(zeta, omega_n, capacity_cost, power_lim_frac, B_p, location, params)

    data_V1 = params.cem_data;

    data_V1.wec_cost( data_V1.wec_cost == 5000 ) = 15000; % fixme: this is just placeholder for testing

    %{
    CEM_co2 = interpn(zeta_data, omega_n_data, cap_cost_data, power_lim_data, CO2_data,...
        zeta,omega_n,capacity_cost,power_limit);
    CEM_capacity = interpn(zeta_data, omega_n_data, cap_cost_data, power_lim_data, capacity_data,...
        zeta,omega_n,capacity_cost,power_limit);
    %}

    
    %if zeta == 0.05 && omega_n == 0.4 && strcmp(location,'ISONE')

    CEM_data = findNearestRow_interp(zeta, omega_n, capacity_cost, power_lim_frac*1000, data_V1);

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
            capacity_cost_pct_incr = (capacity_cost - 725) / 725;
    
            %CEM_CO2 =  7.551749e6 * (1 + capacity_cost_pct_incr * co2_slope);
            %CEM_wec_capacity = 10.201e3 * (1 - capacity_cost_pct_incr * cap_slope);
            %CEM_grid_cost = 2.468040544e9 * (1 - capacity_cost_pct_incr * cost_slope);

            CEM_CO2 = CEM_data.carbon;
            CEM_wec_capacity = CEM_data.wave_capacity;
            CEM_grid_cost = CEM_data.system_cost;

            
        else 
            % not in bounds of model

            %cheapest_cost_with_data
            % capacity_cost
            error('WEC is too cheap, no CEM data here.')
        end
    else
        error('CEM results not available for this zeta, w_n, and location')
    end

end