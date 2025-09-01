classdef DescFcns < GenericAnalysis
    %DESCFCNS Analysis class for describing function figures
    %   Generates drag and saturation describing function figures
    
    properties
        fig_names = {'drag_desc_fcn', 'saturation_desc_fcn'};
        tab_names = {};
    end
    
    methods
        
        function intermed_result_struct = analysis_fcn(~)
            % Run describing function demo
            sin_desc_fcn_demo()
            
            % Store figure numbers for post-processing
            intermed_result_struct.drag_fig_number = gcf().Number;
        end
        
        function [fig_array,...
                 tab_array_display,...
                 tab_array_latex,...
                 end_result_struct] = post_process_fcn(intermed_result_struct)
            
            drag_fig_num = intermed_result_struct.drag_fig_number;
            
            fig_array = [figure(drag_fig_num), figure(drag_fig_num - 2)];
            
            tab_array_display = {};
            tab_array_latex = {};
            
            end_result_struct.desc_fcn_analysis_complete = true;
        end
        
    end
end
