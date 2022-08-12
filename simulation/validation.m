% Run this validation by typing runtests('validation') in the command window
clc

[~,~,~,tab] = validate_nominal_RM3()

%% Feasible
[feasible,failed] = validate_nominal_RM3();
failed
assert( feasible )

%% Mass within 10 percent
[~,~,pct_error] = validate_nominal_RM3();
err = [pct_error.mass_f pct_error.mass_vc pct_error.mass_rp pct_error.mass_tot];
assert( all(err < 0.1) )

%% Cost within 5 percent
[~,~,pct_error] = validate_nominal_RM3();
err = [pct_error.capex pct_error.opex];
assert( all(err < 0.05) )

%% Power within 10 percent
[~,~,pct_error] = validate_nominal_RM3();
err = [pct_error.power_avg pct_error.power_max];
assert( all(err < 0.1) )

%% Force within 10 percent
[~,~,pct_error] = validate_nominal_RM3();
err = [pct_error.force_heave pct_error.FOS_b];
assert( all(err < 0.1) )

%% LCOE within 10 percent
[~,~,pct_error] = validate_nominal_RM3();
err = pct_error.LCOE;
assert( all(err < 0.1) )

%%%%%% function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [feasible,failed,pct_error,tab] = validate_nominal_RM3()
    p = parameters();
    p.N_WEC = 1;
    p.power_max = 286000;
    p.LCOE_max = 10; % set large max LCOE to avoid failing feasibility check
    b = var_bounds(p); 
    
    X = [b.X_noms; 1];
    
    [~, ~, B, FOS1Y, FOS2Y, FOS3Y, FOS_buckling, ...
                                GM, P_elec, D_d, ~, g, simulated] = simulation(X,p);
    
    FOS = min([FOS1Y FOS2Y FOS3Y FOS_buckling]);
    [feasible,failed] = is_feasible(B,FOS,GM,P_elec,D_d,g(16),g(17),g(18),p);
    
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
                    [~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, tmp(j)] = simulation(X,p);
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

        if nargout > 3
            % create combined struct
            results = simulated;
            results(2) = actual;
            results(3) = pct_error;
            tab = struct2table(results, 'RowNames',{'Simulation','RM3 actual','Error'});
        end
    end
end
