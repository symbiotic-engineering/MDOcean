% find the value of F_max for the nominal RM3 
% that is just high enough so there is no saturation

function F_max_nom = find_nominal_inputs(b, mode, display_on)

    % setup
    p = parameters(mode);
    p.N_WEC = 1;

    % optimize
    F_max_nom = fmincon(@(F_max)errFunc(F_max,p,b), b.F_max_nom,[],[],[],[],b.F_max_min,b.F_max_max);

    if display_on
        % check feasibility
        X = [b.X_noms; 1];
        idx_F = strcmp(b.var_names,'F_max');
        X(idx_F) = F_max_nom;
        [LCOE, P_var, ~, g] = simulation(X,p)
        [feasible,~,failed] = is_feasible(g, X, p, b)

        % display x output
        array2table(F_max_nom,'VariableNames',{'F_max (1e6 N)'})
        
        % display y output
        [~,y] = errFunc(x,y_desired,p,b);
        results = round([y; y_desired] ./ 1000,3,'significant');
        array2table(results,'RowNames',{'Sim Output','RM3 Actual'},...
            'VariableNames','Max Powertrain Force Error (-)')
    end

end

function err = errFunc(F_max_in,p,b)
    X = [b.X_noms; 1];
    idx_F = strcmp(b.var_names,'F_max');
    X(idx_F) = F_max_in;

    [~, ~, ~, g] = simulation(X, p);
    idx_F_constr = strcmp(b.constraint_names,'irrelevant_max_force');
    err = g(idx_F_constr)^2;
end
