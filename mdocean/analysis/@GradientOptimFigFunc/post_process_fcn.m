function [fig_array,...
                 tab_array_display,...
                 tab_array_latex,...
                 end_result_struct] = post_process_fcn(intermed_result_struct)
    
    p = intermed_result_struct.p;
    b = intermed_result_struct.b;
    which_objs = intermed_result_struct.which_objs;
    Xs_opt = intermed_result_struct.Xs_opt;
    objs_opt = intermed_result_struct.objs_opt;
    flags = intermed_result_struct.flags;
    lambdas = intermed_result_struct.lambdas;
    grads = intermed_result_struct.grads;
    hesses = intermed_result_struct.hesses;

    [fig_array,tab] = SOO_result_plots(Xs_opt,lambdas,grads,hesses,objs_opt,which_objs,p,b);
    
    fig_array(end+1) = intermed_result_struct.convergence_plot;

    tab_array_display = {tab};
    tab_array_latex = {tab};
    
    end_result_struct.convergence_achieved = all(flags > 0);
end
