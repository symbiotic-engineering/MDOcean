function intermed_result_struct = analysis_fcn(~,~)
    % Run MEEM validation analysis
    [figPotMatch, figVelMatch, figSparsity, figHydroCoeff] = validate_MEEM();
    b_vector_fig = b_inf_numeric();
    [fig_convergence_vs_omega,fig_convergence_vs_NMK] = convergence_study();
    
    % Store results for post-processing
    intermed_result_struct.figPotMatch = figPotMatch;
    intermed_result_struct.figVelMatch = figVelMatch;
    intermed_result_struct.figSparsity = figSparsity;
    intermed_result_struct.figHydroCoeff = figHydroCoeff;
    intermed_result_struct.b_vector_fig = b_vector_fig;
    intermed_result_struct.fig_convergence_vs_omega = fig_convergence_vs_omega;
    intermed_result_struct.fig_convergence_vs_NMK = fig_convergence_vs_NMK;
end
