function [net_eco_value] = net_value(CEM_CO2,eco_cost_total,social_cost_carbon,env_cost_per_wec,AEP_matrix,p, force_matrix)

%enviro
CEM_CO2_kg = CEM_CO2 * 1e3; % ton to kg
eco_value_total = CEM_CO2_kg * social_cost_carbon;
net_eco_value = (eco_value_total - eco_cost_total);

%grid cem
levelized_eco_value_matrix = AEP_matrix .* p.marginal_carbon * p.eco_cost_carbon; % kWh/year * kgCO2/kWh * $/kgCO2 = $/year
levelized_eco_cost = env_cost_per_wec * p.N_WEC; % fixme this isn't levelized
levelized_eco_cost_matrix = levelized_eco_cost * force_matrix;

[LCOE_matrix_eco,LCOE_eco,...
 LVOE_matrix_eco,LVOE_eco,...
 NVOE_matrix_eco,NVOE_eco,...
 net_value_matrix_eco,net_value_eco,...
 VCR_matrix_eco, VCR_eco] = get_LXOE_NVOE_VCR(levelized_eco_cost_matrix, ...
                                              levelized_eco_value_matrix, ...
                                              AEP_matrix);
net_eco_value = mysum(net_value_matrix_eco);

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