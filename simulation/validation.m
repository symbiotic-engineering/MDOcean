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
pct_error(1:3)
assert( all(pct_error(1:3) < 0.01) )

%% Cost within 5 percent
[~,~,pct_error] = validate_nominal_RM3();
disp('Cost percent error:')
pct_error(4:5)
assert( all(pct_error(4:5) < 0.05) )

%% Power within 10 percent
[~,~,pct_error] = validate_nominal_RM3();
disp('Power percent error:')
pct_error(6)
assert( pct_error(6) < 0.1 )

%% Force within 10 percent
[~,~,pct_error] = validate_nominal_RM3();
disp('Force percent error:')
pct_error(7:8)
assert( all(pct_error(7:8) < 0.1) )

%% LCOE within 10 percent
[~,~,pct_error] = validate_nominal_RM3();
disp('LCOE percent error:')
pct_error(9)
assert( pct_error(9) < 0.1 )


%%%%%% function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [feasible,failed,pct_error,tab] = validate_nominal_RM3()
    p = parameters();
    p.N_WEC = 1;
    b = var_bounds(p);
    
    
    X = [b.X_noms; 1];
    
    [LCOE, P_var, B, FOS1Y, FOS2Y, FOS3Y, FOS_buckling, ...
                                GM, P_elec, D_d, ~, g, val] = simulation(X,p);
    
    FOS = min([FOS1Y FOS2Y FOS3Y FOS_buckling]);
    [feasible,failed] = is_feasible(B,FOS,GM,P_elec,D_d,g(16),g(17),g(18),p);
    
    if nargout > 2
        mass_f_nom = 208e3;
        mass_vc_nom = 224e3;
        mass_rp_nom = 245e3;
        capex_nom = 17e6;
        opex_nom = 1.2e6;
        power_nom = 286e3;
        F_heave_nom = 8500e3;
        FOS_b_nom = 3;
        LCOE_nom = 0.75;
        
        actual = [mass_f_nom mass_vc_nom mass_rp_nom capex_nom opex_nom ...
                power_nom F_heave_nom FOS_b_nom LCOE_nom];
        simulated = [val min(FOS_buckling) LCOE];
        
        pct_error = abs(simulated-actual) ./ actual;

        if nargout > 3
            results = [actual; simulated; pct_error];
            tab = array2table(results, 'RowNames',{'RM3 actual','Simulation','Error'},...
                'VariableNames',{'Float mass','VC mass','RP mass',...
                'Capex','Opex','Power','Heave force','Buckling FOS','LCOE'});
        end
    end
end
