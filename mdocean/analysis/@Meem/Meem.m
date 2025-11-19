classdef Meem < GenericAnalysis
    %MEEM Analysis class for MEEM validation figures
    %   Generates MEEM regions, sparsity, validation, matching, and convergence figures

    properties
        fig_names = {'meem_regions', 'meem_sparsity', 'meem_validation', 'meem_matching', 'meem_convergence', 'asymptotic_b_vector'};
        tab_names = {};
    end

    methods (Static)
        intermed_result_struct = analysis_fcn(~,~)

        [fig_array,...
                 tab_array_display,...
                 tab_array_latex,...
                 end_result_struct] = post_process_fcn(intermed_result_struct)
    end
end
