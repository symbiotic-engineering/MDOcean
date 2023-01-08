function [feasible,failed,simulated,actual,tab] = validate_nominal_RM3()
    p = parameters();
    p.N_WEC = 1;
    p.power_max = 286000;
    p.LCOE_max = 10; % set large max LCOE to avoid failing feasibility check
    b = var_bounds(p); 
    
    X = [b.X_noms; 1];
    
    [~, ~, ~, g, simulated] = simulation(X,p);
    
    [feasible,failed] = is_feasible(g, b);

    if nargout > 2
        actual = validation_inputs();
        tiledlayout(1,3)
        fields = fieldnames(actual);
        for i = 1:length(fields)
            if any(strcmp(fields{i},{'capex','opex','LCOE'}))
                % for economic validation, sweep N_WEC
                N_WEC = [1 10 50 100];
                tmp = simulated;
                for j = 2:length(N_WEC)
                    p.N_WEC = N_WEC(j); 
                    [~, ~, ~, ~, tmp(j)] = simulation(X,p);
                end
                simulated.(fields{i}) = [tmp.(fields{i})];  
                
                nexttile
                semilogx(N_WEC,simulated.(fields{i}),N_WEC,actual.(fields{i}))
                xlabel('N_{WEC}')
                title((fields{i}))
                legend('Simulated','Actual')
            end
            sim = simulated.(fields{i}); 
            act = actual.(fields{i});
            pct_error.(fields{i}) = abs(sim-act) ./ act;
        end
        improvePlot

        if nargout > 4
            % create combined struct
            results = rmfield(simulated, 'power_unsat');
            results(2) = actual;
            results(3) = pct_error;
            tab = struct2table(results, 'RowNames',{'Simulation','RM3 actual','Error'});
        end
   end
end
