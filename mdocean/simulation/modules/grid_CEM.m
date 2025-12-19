function [CEM_CO2,CEM_grid_cost,AEP_matrix,force_matrix] = grid_CEM(F_heave_mat,capex,opex,P_matrix_elec,sim_LCOE)

p = parameters();
in = p;

%% power/energy
% power
P_weighted = P_matrix_elec .* p.JPD / 100 * in.eff_array; % taken from dynamics.m and LCOE_from_capex_design_power.m
P_avg = mysum(P_weighted);

% annual energy
hr_per_yr = 8766;
AEP_matrix = P_weighted * in.N_WEC * hr_per_yr / 1000; % W to kWh per year, all wecs
AEP = P_avg * in.N_WEC * hr_per_yr / 1000; % kWh per year

%% force for scaling
%scale the cost by the amount of steel needed at a given sea state
%capex equals the sum of the entries of capex_matrix
factor = 1./mysum(F_heave_mat);
force_matrix = F_heave_mat.*factor;

%% econ
levelized_cost = in.FCR * capex + opex; % $/year
levelized_cost_matrix_unscaled = in.FCR * capex * force_matrix + opex; % $/year
cost_matrix_scale_factor = levelized_cost / mysum(levelized_cost_matrix_unscaled);
levelized_cost_matrix = levelized_cost_matrix_unscaled * cost_matrix_scale_factor;

levelized_value_matrix = AEP_matrix .* p.marginal_price; % kWh per year * $/kWh = $/year

[LCOE_matrix,LCOE,...
 LVOE_matrix,LVOE,...
 NVOE_matrix,NVOE,...
 net_value_matrix,net_value,...
 VCR_matrix, VCR] = get_LXOE_NVOE_VCR(levelized_cost_matrix, ...
                                      levelized_value_matrix, ...
                                      AEP_matrix);

% lcoe check
assert(ismembertol(sim_LCOE,LCOE));

CEM_CO2 = in.marginal_carbon;
CEM_grid_cost = net_value;

%% wave conditions
[T_mdocean,Hs_mdocean] = meshgrid(p.T,p.Hs);



function [LCOE_matrix,LCOE,...
          LVOE_matrix,LVOE,...
          NVOE_matrix,NVOE,...
          net_value_matrix,net_value,...
          VCR_matrix, VCR] = get_LXOE_NVOE_VCR(levelized_cost_matrix, ...
                                               levelized_value_matrix, ...
                                               AEP_matrix)
    
    % lcoe = levelized cost / annual energy production
    [LCOE_matrix, LCOE] = LXOE(levelized_cost_matrix,  AEP_matrix);

    % lvoe = levelized value / annual energy production
    [LVOE_matrix, LVOE] = LXOE(levelized_value_matrix, AEP_matrix);

    % net value, NVOE, and value cost ratio matrices
    net_value_matrix = levelized_value_matrix - levelized_cost_matrix;
    NVOE_matrix      = LVOE_matrix            -  LCOE_matrix;
    VCR_matrix       = LVOE_matrix           ./  LCOE_matrix;

    % net value, NVOE, and value cost ratio scalars
    net_value = mysum(net_value_matrix);
    NVOE = LVOE - LCOE;
    VCR  = LVOE / LCOE;
    
end

function [LXOE_matrix,LXOE] = LXOE(levelized_X_matrix,AEP_matrix)
    % lxoe = levelized X (where X = cost or value) / annual energy production
    LXOE_matrix = levelized_X_matrix ./ AEP_matrix; % from LCOE_from_capex_design_power.m: $/year / kWh/year = $/kWh
    levelized_X = mysum(levelized_X_matrix);
    AEP = mysum(AEP_matrix);
    LXOE = levelized_X ./ AEP;
end

function out = mysum(in)
    out = sum(in,'all','omitnan');
end



end