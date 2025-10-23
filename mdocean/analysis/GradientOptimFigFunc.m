classdef GradientOptimFigFunc < GenericAnalysis
    %GRADIENTOPTIMFIGFUNC Analysis class for gradient optimization figures
    %   Generates single objective optimization convergence and visualization figures
    
    properties
        fig_names = {'dJ_dx_gradient', 'lagrange_multipliers', 'single_obj_opt_geometry', ...
                    'single_obj_opt_power_matrix', 'single_obj_convergence', 'normalized_gradient','delta_x'};
        tab_names = {};
    end
    
    methods (Static)
        
        function intermed_result_struct = analysis_fcn(p,b)
            % Run gradient optimization with plotting
            [Xs_opt, objs_opt, flags, probs, ...
             lambdas, grads, hesses, vals, figs] = gradient_optim(b.X_start_struct,p,b,1,...
                                                               {@optimplotfval, @(x,~,~)optim_geomviz(x,p,b)},true)

            intermed_result_struct.Xs_opt = Xs_opt;
            intermed_result_struct.objs_opt = objs_opt;
            intermed_result_struct.flags = flags;
            intermed_result_struct.probs = probs;
            intermed_result_struct.lambdas = lambdas;
            intermed_result_struct.grads = grads;
            intermed_result_struct.hesses = hesses;
            intermed_result_struct.vals = vals;
            intermed_result_struct.figs = figs;
        end
        
        function [fig_array,...
                 tab_array_display,...
                 tab_array_latex,...
                 end_result_struct] = post_process_fcn(intermed_result_struct)
            
            % Get the figures created by gradient_optim
            fig_array = intermed_result_struct.figs;
            
            tab_array_display = {};
            tab_array_latex = {};
            
            end_result_struct.convergence_achieved = all(intermed_result_struct.flags > 0);
        end
        
    end
end
