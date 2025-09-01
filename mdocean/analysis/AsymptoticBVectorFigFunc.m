classdef AsymptoticBVectorFigFunc < GenericAnalysis
    %ASYMPTOTICBVECTORFIGFUNC Analysis class for asymptotic b vector figures
    %   Generates asymptotic b vector figure
    
    properties
        fig_names = {'asymptotic_b_vector'};
        tab_names = {};
    end
    
    methods
        
        function intermed_result_struct = analysis_fcn(~)
            % Generate asymptotic b vector figure
            fig = b_inf_numeric();
            
            % Store figure for post-processing
            intermed_result_struct.figure_handle = fig;
        end
        
        function [fig_array,...
                 tab_array_display,...
                 tab_array_latex,...
                 end_result_struct] = post_process_fcn(intermed_result_struct)
            
            fig_array = intermed_result_struct.figure_handle;
            
            tab_array_display = {};
            tab_array_latex = {};
            
            end_result_struct.asymptotic_b_vector_complete = true;
        end
        
    end
end
