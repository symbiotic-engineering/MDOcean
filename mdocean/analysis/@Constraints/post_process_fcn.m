function [fig_array,...
                 tab_array_display,...
                 tab_array_latex,...
                 end_result_struct,...
                 tab_firstrows,...
                 tab_colspecs] = post_process_fcn(intermed_result_struct)
            
            tab = intermed_result_struct.constraint_names;
            
            fig_array = [];
            
            tab_array_display = {tab};
            tab_array_latex = {tab};
            
            tab_firstrows = {[]};
            tab_colspecs = {[]};
            
            end_result_struct.constraints_table = tab;
end
