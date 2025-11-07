function [figs,tab] = SOO_result_plots(Xs_opt,lambdas,grads,hesses,objs_opt,which_objs,p,b)

    num_objectives_to_run = length(which_objs);
    figs_per_obj = 4;
    figs = gobjects(1,num_objectives_to_run*figs_per_obj);

    for i = 1:num_objectives_to_run
        X_opt = Xs_opt(:,i);
        lambda = lambdas(i);
        grad = grads(:,i);
        hess = hesses(:,:,i);
        obj_opt = objs_opt(i);
        which_obj = which_objs(i);
        start_fig_num = (i-1)*figs_per_obj;

        figs(start_fig_num + 1) = plot_power_matrix(X_opt,p,b,b.filename_uuid);
        figs(start_fig_num + 2) = visualize_geometry(X_opt,p);
        figs(start_fig_num + 3) = lagrange_multiplier_bar_chart(b,lambda);
        figs(start_fig_num + 4) = delta_x(X_opt,grad,hess,obj_opt,p,b,which_obj);
    end
    table_data = [Xs_opt(1:end-1,:), b.X_mins, b.X_maxs b.X_noms];
    tab = array2table(table_data,'RowNames',b.var_names(1:end-1),...
            'VariableNames',[strcat("Min ", b.obj_names(which_objs)), {'Min bound','Max bound','Nom'}]);
end
