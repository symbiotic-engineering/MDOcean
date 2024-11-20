function [feasible,failed,simulated,actual,tab] = validate_nominal_RM3(mode)
    p = parameters(mode);
    p.N_WEC = 1;
    p.power_max = 286000;
    p.LCOE_max = 10; % set large max LCOE to avoid failing feasibility check
    b = var_bounds(mode); 
    
    X = [b.X_noms; 1];
    
    [~, ~, ~, g, simulated] = simulation(X,p);
    
    % whether nominal violates constraints
    [feasible,failed] = is_feasible(g, X, p, b);

    % comparison of simulated and actual values
    if nargout > 2
        actual = validation_inputs(mode);
        tiledlayout(1,3)
        fields = fieldnames(actual);

        % for economic validation, sweep N_WEC
        N_WEC = [1 10 50 100];
        simulated_diff_N_WEC = simulated;
        for j = 2:length(N_WEC)
            p.N_WEC = N_WEC(j); 
            [~, ~, ~, ~, simulated_diff_N_WEC(j)] = simulation(X,p);
        end

        % sweep through each field (each item to validate)
        for i = 1:length(fields)
            field = fields{i};

            % if the field is economic, plot vs N_WEC
            if any(strcmp(field,{'capex','opex','LCOE'}))
                simulated.(field) = [simulated_diff_N_WEC.(field)];  
                
                nexttile
                semilogx(N_WEC,simulated.(field),N_WEC,actual.(field))
                xlabel('N_{WEC}')
                title((field))
                legend('Simulated','Actual')
            end

            % if the field is dynamic, plot vs sea state matrix
            if any(strcmp(field,{'P','X','B_p'}))
                % fixme: add plots from power_matrix_compare here
            end

            sim = simulated.(field); 
            act = actual.(field);
            pct_error.(field) = abs(sim-act) ./ act;
        end
        improvePlot

        % simulated and actual in table form 
        if nargout > 4
            % create combined struct
            extra_fields = setdiff(fieldnames(simulated),fieldnames(actual));
            results = rmfield(simulated, extra_fields); % remove fields where there is no actual data
            results(2) = actual;
            results(3) = pct_error;
            tab = struct2table(results, 'RowNames',{'Simulation','RM3 actual','Error'});
        end
   end
end
