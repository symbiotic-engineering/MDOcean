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
