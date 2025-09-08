% find the value of F_max for the nominal RM3 
% that is just high enough so there is no saturation

function F_max_nom = find_nominal_inputs(b, p)

    % optimize: max F_max s.t. g >= 0 (force limit active)
    obj = @(F_max)deal(-F_max, -1);
    const = @(F_max)constrFunc(F_max,p,b);
    x0 = b.F_max_nom;
    A_ineq = []; B_ineq = []; A_eq = []; B_eq = [];
    lb = b.F_max_min;
    ub = b.F_max_max;
    options = optimoptions("fmincon",'SpecifyObjectiveGradient',true,'Display','notify');
    F_max_nom = fmincon(obj, x0, A_ineq, B_ineq, A_eq, B_eq, lb, ub, const, options);

end

function [g_out, g_eq] = constrFunc(F_max_in,p,b)
    X = [b.X_noms; 1];
    idx_F = strcmp(b.var_names,'F_max');
    X(idx_F) = F_max_in;

    [~, ~, g] = simulation(X, p);
    idx_F_constr = strcmp(b.constraint_names,'irrelevant_max_force');
    g_out = -g(idx_F_constr);
    g_eq = [];
end
