function [fig_array,...
                 tab_array_display,...
                 tab_array_latex,...
                 end_result_struct] = post_process_fcn(intermed_result_struct)
            
            tab = array2table([intermed_result_struct.X_mins intermed_result_struct.X_noms intermed_result_struct.X_maxs], ...
                    'VariableNames',{'Mins','Noms','Maxs'}, 'RowNames', intermed_result_struct.var_names(1:end-1));
            
            fig_array = [];
            
            tab_array_display = {tab};
            tab_array_latex = {tab};
            
            end_result_struct.design_vars_table = tab;
end
