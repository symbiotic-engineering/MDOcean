classdef GradientOptimFigFunc < GenericAnalysis
    %GRADIENTOPTIMFIGFUNC Analysis class for gradient optimization figures
    %   Generates single objective optimization convergence and visualization figures
    
    properties
        fig_names = {'dJ_dx_gradient', 'lagrange_multipliers', 'single_obj_opt_geometry', ...
                    'single_obj_opt_power_matrix', 'single_obj_convergence', 'normalized_gradient'};
        tab_names = {};
    end
    
    methods (Static)
        
        function intermed_result_struct = analysis_fcn(p,b)
            % Run gradient optimization with plotting
            gradient_optim(b.X_start_struct,p,b,1,{@optimplotfval, @(x,~,~)optim_geomviz(x,p,b)},true)

            % Store figure numbers for post-processing
            n = gcf().Number;
            intermed_result_struct.final_figure_number = n;
        end
        
        function [fig_array,...
                 tab_array_display,...
                 tab_array_latex,...
                 end_result_struct] = post_process_fcn(intermed_result_struct)
            
            n = intermed_result_struct.final_figure_number;
            
            % Get the figures created by gradient_optim
            fig_array = [figure(n-4), figure(n-3), figure(n-2), figure(n-1), figure(n), figure(n+1)];
            
            tab_array_display = {};
            tab_array_latex = {};
            
            end_result_struct.convergence_achieved = true;
        end
        
    end
end
