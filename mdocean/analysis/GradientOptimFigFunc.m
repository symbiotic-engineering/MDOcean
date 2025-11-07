classdef GradientOptimFigFunc < GenericAnalysis
    %GRADIENTOPTIMFIGFUNC Analysis class for gradient optimization figures
    %   Generates single objective optimization convergence and visualization figures

    properties
        fig_names = {'single_obj_opt_power_matrix',...
                    'single_obj_opt_geometry', ...
                    'lagrange_multipliers',...
                    'delta_x', ...
                    'single_obj_convergence'};
        tab_names = {'single_obj_optim_results'};
    end

    methods (Static)

        function intermed_result_struct = analysis_fcn(p,b)
            which_objs = 1;

            % Run gradient optimization
            [Xs_opt, objs_opt, flags, probs, ...
             lambdas, grads, hesses, vals] = gradient_optim(b.X_start_struct,p,b,which_objs,...
                                                               {@optimplotfval, @(x,~,~)optim_geomviz(x,p,b)},true);

            intermed_result_struct.p = p;
            intermed_result_struct.b = b;
            intermed_result_struct.which_objs = which_objs;
            intermed_result_struct.Xs_opt = Xs_opt;
            intermed_result_struct.objs_opt = objs_opt;
            intermed_result_struct.flags = flags;
            intermed_result_struct.probs = probs;
            intermed_result_struct.lambdas = lambdas;
            intermed_result_struct.grads = grads;
            intermed_result_struct.hesses = hesses;
            intermed_result_struct.vals = vals;
            intermed_result_struct.convergence_plot = gcf();

        end

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

    end
end
