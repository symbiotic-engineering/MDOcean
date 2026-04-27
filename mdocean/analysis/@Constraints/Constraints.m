classdef Constraints < GenericAnalysis
    %CONSTRAINTS Analysis class for constraints table
    %   Generates constraints information table

    properties
        fig_names = {};
        tab_names = {'constraints'};
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
