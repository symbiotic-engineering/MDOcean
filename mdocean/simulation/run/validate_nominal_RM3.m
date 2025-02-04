function [feasible,failed,simulated,actual,tab,fig] = validate_nominal_RM3(mode)
    p = parameters(mode);
    p.N_WEC = 1;
    p.LCOE_max = 10; % set large max LCOE to avoid failing feasibility check
    p.control_type = 'damping';
    b = var_bounds(mode); 
    
    X = [b.X_noms; 1];
    
    [~, ~, ~, g, simulated] = simulation(X,p);
    
    % whether nominal violates constraints
    idx_ignore = strcmp(b.constraint_names,'irrelevant_max_force');
    [feasible,failed] = is_feasible(g, X, p, b, idx_ignore);

    % comparison of simulated and actual values
    if nargout > 2
        actual = validation_inputs(mode);
        fig = figure;
        econ_fields = {'capex','opex','LCOE','capex_design','capex_struct','capex_PTO','J_capex_design'};
        t = tiledlayout(fig,1,length(econ_fields)-1);
        t.TileSpacing = 'compact';
        t.Padding = 'compact';
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

            % if the field is economic (execpt J_capex_design since it's a duplicate), plot vs N_WEC
            if any(strcmp(field,econ_fields)) && ~strcmp(field,'J_capex_design')
                simulated.(field) = [simulated_diff_N_WEC.(field)];  
                
                ax = nexttile(t);
                semilogx(ax, N_WEC,simulated.(field),N_WEC,actual.(field),'--')
                xlabel(ax,'N_{WEC}')
                title(ax,remove_underscores({field}))
                if length(t.Children) == t.GridSize(2)
                    legend(ax,'Simulated','Actual','Location','bestoutside') % only add legend to last plot
                end
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
        set(fig,"Position",[3.4 201 1529.2 600])
        set(ax.Legend,'Position',[0.2671 0.1913 0.0973 0.0902])

        % simulated and actual in table form 
        if nargout > 4
            % create combined struct
            extra_fields = setdiff(fieldnames(simulated),fieldnames(actual));
            results = rmfield(simulated, extra_fields); % remove fields where there is no actual data
            results(2) = actual;
            results(3) = pct_error;
            tab = struct2table(results, 'RowNames',strcat({'Simulation','RM3 actual','Error'},{[' ' mode]}));
        end
   end
end
