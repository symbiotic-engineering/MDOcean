classdef DesignVars < GenericAnalysis
    %DESIGNVARS Analysis class for design variables table
    %   Generates design variables bounds table

    properties
        fig_names = {};
        tab_names = {'design_vars'};
    end

    methods (Static)
        intermed_result_struct = analysis_fcn(~,b)

        [fig_array,...
                 tab_array_display,...
                 tab_array_latex,...
                 end_result_struct] = post_process_fcn(intermed_result_struct)
    end
end
