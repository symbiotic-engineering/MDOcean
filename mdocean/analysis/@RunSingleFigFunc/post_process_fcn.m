function [fig_array,...
                 tab_array_display,...
                 tab_array_latex,...
                 end_result_struct] = post_process_fcn(intermed_result_struct)

    fig_array = intermed_result_struct.created_figs;
    
    tab_array_display = {};
    tab_array_latex = {};
    
    end_result_struct.single_run_complete = true;
end
