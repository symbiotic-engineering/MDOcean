
function net_eco_value = enviro(steel_mass, distance_from_shore, fiberglass_area, ...
                                CEM_CO2, CEM_WEC_capacity, rated_power_per_wec, ...
                                eco_cost_steel, eco_cost_fiberglass, ...
                                eco_cost_distance, social_cost_carbon)

    % units: [steel_kg, distance_miles, fiberglass_meters_sq, CEM_CO2_tonne]

    CEM_CO2_kg = CEM_CO2 * 1e3; % ton to kg
    num_wecs_CEM = CEM_WEC_capacity / rated_power_per_wec;

    % cost and value
    eco_value_total = CEM_CO2_kg * social_cost_carbon;
    eco_cost_per_wec = eco_cost_steel * steel_mass + eco_cost_fiberglass * fiberglass_area + eco_cost_distance * distance_from_shore;
    eco_cost_total = eco_cost_per_wec * num_wecs_CEM;

    % convert eco
    net_eco_value = (eco_value_total - eco_cost_total);
end