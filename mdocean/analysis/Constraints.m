classdef Constraints < GenericAnalysis
    %CONSTRAINTS Analysis class for constraints table
    %   Generates constraints information table
    
    properties
        fig_names = {};
        tab_names = {'constraints'};
    end
    
    methods
        
        function intermed_result_struct = analysis_fcn(obj)
            % Get constraints information
            % Store results for post-processing
            intermed_result_struct.constraint_names = obj.b.constraint_names;
        end
        
        function [fig_array,...
                 tab_array_display,...
                 tab_array_latex,...
                 end_result_struct] = post_process_fcn(intermed_result_struct)
            
            tab = intermed_result_struct.constraint_names';
            
            fig_array = [];
            
            tab_array_display = {tab};
            tab_array_latex = {tab};
            
            end_result_struct.constraints_table = tab;
        end
        
    end
end
