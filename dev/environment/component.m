function optimize_environmental_value()
    % Constants
    s_points = 0.1764;        % euro/kg (steel)
    f_points = 6.3757;        % euro/m^2 (fiberglass)
    d_points = 60.4413;       % euro/mi (distance)
    SCC = 0.133;              % euro/kg CO2
    euro2USD = 1.1;           % currency conversion

    lb = [0, 0, 0, 0];        % lower bounds
    ub = [5000, 100, 1000, 0.01];  % upper bounds (example limits)

    % initialization
    x0 = [1000, 50, 200, 0.005];  % [steel, distance, fiberglass, CEM_output]

    options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'sqp')
    [x_opt, fval] = fmincon(@(x) -net_eco_value(x, s_points, f_points, ...
                    d_points, SCC, euro2USD), ...
                    x0, [], [], [], [], lb, ub, [], options);
    % resultant vector is structured as [steel, distance, fiberglass, CO2]

end

function value = net_eco_value(x, s_points, f_points, d_points, SCC, euro2USD)
    % x = [steel, distance, fiberglass, CEM_output]
    steel = x(1);
    distance = x(2);
    fiberglass = x(3);
    CEM_output = x(4);

    CEM_output_kg = CEM_output * 1e9;

    % cost and value
    eco_value = CEM_output_kg * SCC;
    eco_cost = s_points * steel + f_points * fiberglass + d_points * distance;

    % convert eco
    value = (eco_value - eco_cost) * euro2USD;
end