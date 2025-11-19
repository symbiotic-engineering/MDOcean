classdef ParetoSweep < GenericAnalysis
    %PARETOSWEEP Analysis class for Pareto sweep figures
    %   Generates sweep number of seeds figure

    properties
        fig_names = {'sweep_num_seeds'};
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
