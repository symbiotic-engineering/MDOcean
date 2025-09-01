classdef JpdMultiplyFigFunc < GenericAnalysis
    %JPDMULTIPLYFIGFUNC Analysis class for JPD multiplication figures
    %   Generates JPD multiplication power matrix figure
    
    properties
        fig_names = {'JPD_multiplication'};
        tab_names = {};
    end
    
    methods
        
        function intermed_result_struct = analysis_fcn(obj)
            % Generate JPD multiplication figure
            X = [obj.b.X_noms; 1];
            plot_power_matrix(X, obj.p, obj.b, obj.b.filename_uuid)
            
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
            
            end_result_struct.jpd_multiplication_complete = true;
        end
        
    end
end
