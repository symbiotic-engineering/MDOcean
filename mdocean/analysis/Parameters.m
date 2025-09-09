classdef Parameters < GenericAnalysis
    %PARAMETERS Analysis class for parameters table
    %   Generates parameters information table
    
    properties
        fig_names = {};
        tab_names = {'parameters'};
    end
    
    methods (Static)
        
        function intermed_result_struct = analysis_fcn(~,~)
            % Get parameters information
            [~, params_table] = parameters();
            
            % Store results for post-processing
            intermed_result_struct.params_table = params_table;
        end
        
        function [fig_array,...
                 tab_array_display,...
                 tab_array_latex,...
                 end_result_struct] = post_process_fcn(intermed_result_struct)
            
            tab = intermed_result_struct.params_table;
            
            fig_array = [];
            
            tab_array_display = {tab};
            tab_array_latex = {tab};
            
            end_result_struct.parameters_table = tab;
        end
        
    end
end
