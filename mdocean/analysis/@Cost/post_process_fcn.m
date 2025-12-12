function [fig_array,...
                 tab_array_display,...
                 tab_array_latex,...
                 end_result_struct,...
                 tab_firstrows,...
                 tab_colspecs] = post_process_fcn(intermed_result_struct)
            
            fig_array = [];
            
            tab_array_display = {intermed_result_struct.cost_table};
            tab_array_latex = {intermed_result_struct.cost_table};
            
            tab_firstrows = {[]};
            tab_colspecs = {[]};
            
            end_result_struct.cost_parameters = intermed_result_struct.cost_table;
end
