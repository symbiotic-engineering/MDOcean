classdef DesignSpaceExploration < GenericAnalysis
    %DESIGNSPACEEXPLORATION Analysis class for design space exploration figures
    %   Generates design space exploration experimental figures
    
    properties
        fig_names = {'experiments'};
        tab_names = {};
    end
    
    methods
        
        function intermed_result_struct = analysis_fcn(obj)
            % Run design space exploration experiments
            experiments(obj.p,obj.b)
            
            % Store figure for post-processing
            intermed_result_struct.figure_handle = gcf();
        end
        
        function [fig_array,...
                 tab_array_display,...
                 tab_array_latex,...
                 end_result_struct] = post_process_fcn(intermed_result_struct)
            
            fig_array = intermed_result_struct.figure_handle;
            
            tab_array_display = {};
            tab_array_latex = {};
            
            end_result_struct.design_space_exploration_complete = true;
        end
        
    end
end
