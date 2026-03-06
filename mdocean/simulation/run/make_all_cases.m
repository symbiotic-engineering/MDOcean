function [case_cell,case_desc_cell] = make_all_cases(runOnlyFewSeaStates)

    wecsim_X = [var_bounds('wecsim').X_noms; 1];
    report_X = [var_bounds('report').X_noms; 1];

    wecsim_mb_p_array = make_case_group('wecsim',true,runOnlyFewSeaStates);
    wecsim_sb_p_array = make_case_group('wecsim',false,runOnlyFewSeaStates);
    report_mb_p_array = make_case_group('report',true,runOnlyFewSeaStates);

    case_cell = {wecsim_X, wecsim_mb_p_array;
                wecsim_X, wecsim_sb_p_array;
                report_X, report_mb_p_array};

    case_group_desc = {'drag_off_meem_off','drag_on_meem_off','drag_off_meem_on','drag_on_meem_on'};
    case_desc_cell  = {strcat('wecsim_geom_wecsim_multibody_true_', case_group_desc),...
                       strcat('wecsim_geom_wecsim_multibody_false_', case_group_desc),...
                       strcat('wecsim_geom_report_multibody_true_', case_group_desc),...
                      };

end

function p_array = make_case_group(geometry, multibody, runOnlyFewSeaStates)
    p = parameters(geometry);
    p.use_force_sat = false;
    p.use_power_sat = false;
    p.use_multibody = multibody;

    if runOnlyFewSeaStates
        p.Hs = p.Hs(3:4);
        p.T = p.T(4:5);
        p.JPD = p.JPD(3:4,4:5);
    end

    p_array = make_drag_hydro_param_combos(p);
end

% bonus that I have't checked yet: sweep drag, irregular waves, force saturation

function p_array = make_drag_hydro_param_combos(p_baseline)
    % drag off, wamit coeffs
    p_drag_off_wamit = p_baseline;
    p_drag_off_wamit.C_d_float = 0;
    p_drag_off_wamit.C_d_spar  = 0;
    p_drag_off_wamit.use_MEEM = false;
    
    % drag on, wamit coeffs
    p_drag_on_wamit = p_baseline;
    assert(p_drag_on_wamit.C_d_float>0 && p_drag_on_wamit.C_d_spar>0)
    p_drag_on_wamit.use_MEEM = false;

    % drag off, meem coeffs
    p_drag_off_meem = p_baseline;
    p_drag_off_meem.C_d_float = 0;
    p_drag_off_meem.C_d_spar  = 0;
    p_drag_off_meem.use_MEEM = true;

    % drag on, meem coeffs
    p_drag_on_meem = p_baseline;
    assert(p_drag_on_meem.C_d_float>0 && p_drag_on_meem.C_d_spar>0)
    p_drag_on_meem.use_MEEM = true;
    
    % array of all four
    p_array = [p_drag_off_wamit
               p_drag_on_wamit
               p_drag_off_meem
               p_drag_on_meem];

end