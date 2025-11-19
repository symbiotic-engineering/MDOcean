classdef Cost < GenericAnalysis
    %COST Analysis class for cost analysis tables
    %   Generates cost parameters table

    properties
        fig_names = {};
        tab_names = {'cost_parameters'};
    end

    methods (Static)
        intermed_result_struct = analysis_fcn(~,~)

        [fig_array,...
                 tab_array_display,...
                 tab_array_latex,...
                 end_result_struct] = post_process_fcn(intermed_result_struct)
    end
end
