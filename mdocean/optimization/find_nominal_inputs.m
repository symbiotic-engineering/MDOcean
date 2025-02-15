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
        X(5) = F_max_nom;
        [LCOE, P_var, ~, g] = simulation(X,p)
        [feasible,~,failed] = is_feasible(g, X, p, b)

        % display x output
        array2table(F_max_nom,'VariableNames',{'F_max (1e6 N)','B_p (1e6 Ns/m)','w_n (rad/s)'})
        
        % display y output
        [~,y] = errFunc(x,y_desired,p,b);
        results = round([y; y_desired] ./ [1000 1000 1],3,'significant');
        array2table(results,'RowNames',{'Sim Output','RM3 Actual'},...
            'VariableNames','Max Powertrain Force Error (-)')
    end

end

function err = errFunc(F_max_in,p,b)
    X = [b.X_noms; 1];
    X(5) = F_max_in;
    [~, ~, ~, g] = simulation(X, p);
    err = g(12)^2;
end
