classdef Comparison < GenericAnalysis
    %COMPARISON Analysis class for comparison figures and tables
    %   Generates overlaid geometry, hydro coeffs, probability CDF figures and optimal design tables

    properties
        fig_names = {'overlaid_geometry', 'overlaid_hydro_coeffs', 'probability_CDF', 'comparison_power_matrix'};
        tab_names = {'optimal_design_vars', 'optimal_outputs'};
    end

    methods (Static)
        intermed_result_struct = analysis_fcn(p,b)

        [fig_array,...
         tab_array_display,...
         tab_array_latex,...
         end_result_struct,...
         tab_firstrows,...
         tab_colspecs] = post_process_fcn(intermed_result_struct)
    end
end
