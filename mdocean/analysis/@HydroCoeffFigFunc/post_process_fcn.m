function [fig_array,...
                 tab_array_display,...
                 tab_array_latex,...
                 end_result_struct,...
                 tab_firstrows,...
                 tab_colspecs] = post_process_fcn(intermed_result_struct)
            
            % Return created figures (validation happens in GenericAnalysis)
            fig_array = intermed_result_struct.created_figs;
            
            tab_array_display = {};
            tab_array_latex = {};
            
            tab_firstrows = {};
            tab_colspecs = {};
            
            end_result_struct.hydro_coeff_analysis_complete = true;
end
