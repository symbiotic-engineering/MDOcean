% Run this validation by typing runtests('validation') in the command window
clc

%% Feasible
[feasible,failed] = validate_nominal_RM3();
failed
assert( feasible )

%% Mass within 1 percent
[~,~,pct_error] = validate_nominal_RM3();
pct_error(1:3)
assert( all(pct_error(1:3) < 0.01) )

%% Cost within 5 percent
[~,~,pct_error] = validate_nominal_RM3();
pct_error(4:5)
assert( all(pct_error(4:5) < 0.05) )

%% Power within 10 percent
[~,~,pct_error] = validate_nominal_RM3();
pct_error(6)
assert( pct_error(6) < 0.1 )

%% Force within 10 percent
[~,~,pct_error] = validate_nominal_RM3();
pct_error(7:9)
assert( all(pct_error(7:9) < 0.1) )

%% LCOE within 10 percent
[~,~,pct_error] = validate_nominal_RM3();
pct_error(10)
assert( pct_error(10) < 0.1 )


%%%%%% function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [feasible,failed,pct_error] = validate_nominal_RM3()
    p = parameters();
    p.N_WEC = 1;
    b = var_bounds(p);
    
    
    X = [b.X_noms; 1];
    
    [LCOE, P_var, B, FOS1Y, FOS2Y, FOS3Y, FOS_buckling, ...
                                GM, P_elec, D_d, ~, g, val] = simulation(X,p);
    
    FOS = min([FOS1Y FOS2Y FOS3Y FOS_buckling]);
    [feasible,failed] = is_feasible(B,FOS,GM,P_elec,D_d,g(16),g(17),g(18),p);
    
    if nargout > 1
        mass_f_nom = 208e3;
        mass_vc_nom = 224e3;
        mass_rp_nom = 245e3;
        capex_nom = 17e6;
        opex_nom = 1.2e6;
        power_nom = 286e3;
        F_heave_nom = 8500e3;
        F_ptrain_nom = 8500e3; % fixme
        FOS_b_nom = 3;
        LCOE_nom = 0.75;
        
        actual = [mass_f_nom mass_vc_nom mass_rp_nom capex_nom opex_nom ...
                power_nom F_heave_nom F_ptrain_nom FOS_b_nom LCOE_nom];
        simulated = [val min(FOS_buckling) LCOE];
        
        pct_error = abs(simulated-actual) ./ actual;
    end
end
