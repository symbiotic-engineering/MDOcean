classdef DesignSpaceExploration < GenericAnalysis
    %DESIGNSPACEEXPLORATION Analysis class for design space exploration figures
    %   Generates design space exploration experimental figures
    
    properties
        fig_names = {'experiments'};
        tab_names = {};
    end
    
    methods (Static)
        
        function intermed_result_struct = analysis_fcn(p,b)
            % Run design space exploration experiments
            figs = experiments(p,b);

            % Store figure for post-processing
            intermed_result_struct.figs = figs;
        end
        
        function [fig_array,...
                 tab_array_display,...
                 tab_array_latex,...
                 end_result_struct] = post_process_fcn(intermed_result_struct)
            
            fig_array = intermed_result_struct.figs;
            
            tab_array_display = {};
            tab_array_latex = {};
            
            end_result_struct.design_space_exploration_complete = true;
        end
        
    end
end
