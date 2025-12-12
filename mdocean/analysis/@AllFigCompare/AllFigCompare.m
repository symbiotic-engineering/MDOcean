classdef AllFigCompare < GenericAnalysis
    %ALLFIGCOMPARE Analysis class for runtime comparison bar chart
    %   Generates runtime bar chart comparing all figures

    properties
        fig_names = {'runtime_bar_chart'};
        tab_names = {};
    end

    methods (Static)
        intermed_result_struct = analysis_fcn(~,~)

        [fig_array,...
                 tab_array_display,...
                 tab_array_latex,...
                 end_result_struct,...
                 tab_firstrows,...
                 tab_colspecs] = post_process_fcn(intermed_result_struct)
    end
end
