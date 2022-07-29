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
err = [pct_error.capex pct_error.opex];
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
% [~,~,pct_error] = validate_nominal_RM3();
% disp('LCOE percent error:')
% err = pct_error.LCOE;
% assert( err < 0.1 )
validate_econ_scaling()

%%%%%% function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [feasible,failed,pct_error,tab] = validate_nominal_RM3()
    p = parameters();
    p.N_WEC = 1;
    p.power_max = 286000;
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

function [] = validate_econ_scaling()
    p = parameters();
    p.power_max = 286000;
    b = var_bounds(p);
    X = [b.X_noms; 1];
    
    N_WEC = [1 10 50 100];
    for i = 1:length(N_WEC)
        p.N_WEC = N_WEC(i); 
        [~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, out_struct] = simulation(X,p);
        if i==1
            simulated = out_struct;
        else
            simulated(i) = out_struct;
        end
    end

    actual = validation_inputs();

    figure
    subplot 131
    semilogx(N_WEC,[simulated.LCOE],N_WEC,actual.LCOE)
    xlabel('N_{WEC}')
    title('LCOE')
    legend('Simulated','Actual')

    subplot 132
    semilogx(N_WEC,[simulated.capex],N_WEC,actual.capex)
    xlabel('N_{WEC}')
    title('Capex')
    legend('Simulated','Actual')

    subplot 133
    semilogx(N_WEC,[simulated.opex],N_WEC,actual.opex)
    xlabel('N_{WEC}')
    title('Opex')
    legend('Simulated','Actual')
    improvePlot
end
