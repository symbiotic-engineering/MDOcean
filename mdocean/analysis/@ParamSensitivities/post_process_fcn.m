function [fig_array,...
                 tab_array_display,...
                 tab_array_latex,...
                 end_result_struct,...
                 tab_firstrows,...
                 tab_colspecs] = post_process_fcn(intermed_result_struct)
            
            ratios = intermed_result_struct.ratios;
            LCOE = intermed_result_struct.LCOE;
            LCOE_nom = intermed_result_struct.LCOE_nom;
            P_var = intermed_result_struct.P_var;
            Pvar_nom = intermed_result_struct.Pvar_nom;
            param_names = intermed_result_struct.param_names;
            num_DVs = intermed_result_struct.num_DVs;
            X_LCOE = intermed_result_struct.X_LCOE;
            X_LCOE_nom = intermed_result_struct.X_LCOE_nom;
            dvar_names = intermed_result_struct.dvar_names;
            X_Pvar = intermed_result_struct.X_Pvar;
            X_Pvar_nom = intermed_result_struct.X_Pvar_nom;
            slope_LCOE_norm = intermed_result_struct.slope_LCOE_norm;
            slope_Pvar_norm = intermed_result_struct.slope_Pvar_norm;
            slope_X_LCOE_norm = intermed_result_struct.slope_X_LCOE_norm;
            slope_X_Pvar_norm = intermed_result_struct.slope_X_Pvar_norm;
            param_table = intermed_result_struct.param_table;
            
            par_J_par_p_local = intermed_result_struct.par_J_par_p_local;
            dJ_star_dp_quad_local = intermed_result_struct.dJ_star_dp_quad_local;
            dJ_star_dp_lin_local = intermed_result_struct.dJ_star_dp_lin_local;
            dJstar_dp_global = intermed_result_struct.dJstar_dp_global;
            par_x_star_par_p_local = intermed_result_struct.par_x_star_par_p_local;
            par_x_star_par_p_global = intermed_result_struct.par_x_star_par_p_global;
            delta_p_change_activity_local = intermed_result_struct.delta_p_change_activity_local;
            delta_p_change_activity_global = intermed_result_struct.delta_p_change_activity_global;
            b = intermed_result_struct.b;
            
            figs_global = global_sens_plots(ratios,LCOE,LCOE_nom,P_var,Pvar_nom,param_names,num_DVs,...
                                    X_LCOE,X_LCOE_nom,dvar_names,X_Pvar,X_Pvar_nom,...
                                    slope_LCOE_norm,slope_Pvar_norm,slope_X_LCOE_norm,...
                                    slope_X_Pvar_norm,param_table);

            figs_local_global = all_sens_plots(par_J_par_p_local, dJ_star_dp_quad_local, ...
                                 dJ_star_dp_lin_local, dJstar_dp_global, ...
                                 par_x_star_par_p_local, par_x_star_par_p_global, ...
                                 delta_p_change_activity_local, ...
                                 delta_p_change_activity_global, ...
                                 param_names, dvar_names, b);

            % fig_local_global has 6 and figs_global has 6+num_dvs=18, so 24 total
            fig_array = [figs_local_global,figs_global];
            
            tab_array_display = {param_table};
            tab_array_latex = {param_table};
            
            tab_firstrows = {[]};
            tab_colspecs = {[]};
            
            end_result_struct.sensitivity_analysis_complete = true;
            end_result_struct.runtime_post_optim = intermed_result_struct.runtime_post_optim;
            end_result_struct.runtime_re_optim = intermed_result_struct.runtime_re_optim;
end
