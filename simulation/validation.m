% Run this validation by typing runtests('validation') in the command window
clc

[~,~,~,tab] = validate_nominal_RM3()

%% Feasible
[feasible,failed] = validate_nominal_RM3();
failed
assert( feasible )

%% Mass within 1 percent
[~,~,pct_error] = validate_nominal_RM3();
disp('Mass percent error:')
err = [pct_error.mass_f pct_error.mass_vc pct_error.mass_rp];
assert( all(err < 0.01) )

%% Cost within 5 percent
[~,~,pct_error] = validate_nominal_RM3();
disp('Cost percent error:')
err = [pct_err.capex pct_error.opex];
assert( all(err < 0.05) )

%% Power within 10 percent
[~,~,pct_error] = validate_nominal_RM3();
disp('Power percent error:')
err = [pct_error.power_avg pct_error.power_max];
assert( all(err < 0.1) )

%% Force within 10 percent
[~,~,pct_error] = validate_nominal_RM3();
disp('Force percent error:')
err = [pct_error.force_heave pct_error.FOS_b];
assert( all(err < 0.1) )

%% LCOE within 10 percent
[~,~,pct_error] = validate_nominal_RM3();
disp('LCOE percent error:')
err = pct_error.LCOE;
assert( err < 0.1 )


%%%%%% function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [feasible,failed,pct_error,tab] = validate_nominal_RM3()
    p = parameters();
    p.N_WEC = 1;
    b = var_bounds(p);
    
    
    X = [b.X_noms; 1];
    
    [~, ~, B, FOS1Y, FOS2Y, FOS3Y, FOS_buckling, ...
                                GM, P_elec, D_d, ~, g, simulated] = simulation(X,p);
    
    FOS = min([FOS1Y FOS2Y FOS3Y FOS_buckling]);
    [feasible,failed] = is_feasible(B,FOS,GM,P_elec,D_d,g(16),g(17),g(18),p);
    
    if nargout > 2
        actual = validation_inputs();
        
        fields = fieldnames(actual);
        for i = 1:length(fields)
            sim = simulated.(fields{i});
            act = actual.(fields{i});
            pct_error.(fields{i}) = abs(sim-act) ./ act;
        end
        
        if nargout > 3
            % create combined struct
            results = simulated;
            results(2) = actual;
            results(3) = pct_error;
            tab = struct2table(results, 'RowNames',{'Simulation','RM3 actual','Error'});
        end
    end
end
