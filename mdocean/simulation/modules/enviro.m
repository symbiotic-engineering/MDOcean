function [eco_cost_per_wec,eco_cost_total] = enviro(steel_mass, distance_from_shore, fiberglass_area, ...
                                rated_power_per_wec, ...
                                eco_cost_steel, eco_cost_fiberglass, ...
                                eco_cost_distance, social_cost_carbon)


    % units: [steel_kg, distance_miles, fiberglass_meters_sq, CEM_CO2_tonne]
    num_wecs_CEM = 100; %CEM_WEC_capacity / rated_power_per_wec;

    % cost
    eco_cost_per_wec = eco_cost_steel * steel_mass + eco_cost_fiberglass * fiberglass_area + eco_cost_distance * distance_from_shore;
    eco_cost_total = eco_cost_per_wec * num_wecs_CEM;
end