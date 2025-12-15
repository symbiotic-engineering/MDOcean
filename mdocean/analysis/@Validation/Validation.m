classdef Validation < GenericAnalysis
    %VALIDATION Analysis class for validation figures and tables
    %   Generates cost vs N WEC validation figure and validation table

    properties
        fig_names = {'cost_vs_N_WEC'};
        tab_names = {'validation'};
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
