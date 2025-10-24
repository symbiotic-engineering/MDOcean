classdef LocationSensitivity < GenericAnalysis
    %LOCATIONSENSITIVITY Analysis class for location sensitivity tables
    %   Generates location sensitivity analysis table
    
    properties
        fig_names = {'location_power_matrix', 'location_probability_PDF'};
        tab_names = {'location_sensitivity'};
    end
    
    methods (Static)
        
        function intermed_result_struct = analysis_fcn(p,b)
            % Run location sensitivity analysis
            [tab, h_pm, h_prob] = location_sensitivity(p,b);
            
            % Store results for post-processing
            intermed_result_struct.location_table = tab;
            intermed_result_struct.h_power_matrix = h_pm;
            intermed_result_struct.h_probability_PDF = h_prob;
        end
        
        function [fig_array,...
                 tab_array_display,...
                 tab_array_latex,...
                 end_result_struct] = post_process_fcn(intermed_result_struct)
            
            tab = intermed_result_struct.location_table;
            location_flags = str2double(tab(strcmp(tab.Row,'flag'),:).Variables);
            
            idx_remove = ismember(tab.Row,{'flag','Optimal Material index'});
            tablatex = tab(~idx_remove,:);
            
            fig_array = [intermed_result_struct.h_power_matrix, intermed_result_struct.h_probability_PDF];
            
            tab_array_display = {tab};
            tab_array_latex = {tablatex};
            
            end_result_struct.location_sensitivity_complete = true;
            end_result_struct.location_flags = location_flags;
            end_result_struct.location_table = tab;
        end
        
    end
end
