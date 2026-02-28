classdef DesignSpaceExploration < GenericAnalysis
    %DESIGNSPACEEXPLORATION Analysis class for design space exploration figures
    %   Generates design space exploration experimental figures

    properties
        fig_names = {'experiments_pareto','experiments_ratios'};
        tab_names = {'experiments_results'};
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
