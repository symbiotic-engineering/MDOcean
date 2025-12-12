function [fig_array,...
                 tab_array_display,...
                 tab_array_latex,...
                 end_result_struct,...
                 tab_firstrows,...
                 tab_colspecs] = post_process_fcn(intermed_result_struct)

    fig_array = [intermed_result_struct.fig1, intermed_result_struct.fig2];
    
    tab_array_display = {};
    tab_array_latex = {};
    
    tab_firstrows = {};
    tab_colspecs = {};
    
    end_result_struct.force_saturation_analysis_complete = true;
end
