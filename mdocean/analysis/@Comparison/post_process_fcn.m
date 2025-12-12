function [fig_array,...
                 tab_array_display,...
                 tab_array_latex,...
                 end_result_struct,...
                 tab_firstrows,...
                 tab_colspecs] = post_process_fcn(intermed_result_struct)
                 
    p = intermed_result_struct.p;
    b = intermed_result_struct.b;
    Xs = intermed_result_struct.Xs;
    vals = intermed_result_struct.vals;
    [DV_table, out_table, fig_array] = compare(p,b,Xs,vals);
    
    tab_array_display = {DV_table, ...
                         out_table};
    tab_array_latex = tab_array_display;
    
    tab_firstrows = {[], []};
    tab_colspecs = {[], []};
    
    end_result_struct.comparison_complete = true;
end
