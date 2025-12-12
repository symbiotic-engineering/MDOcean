function [fig_array,...
                 tab_array_display,...
                 tab_array_latex,...
                 end_result_struct,...
                 tab_firstrows,...
                 tab_colspecs] = post_process_fcn(intermed_result_struct)
            
            figs = pareto_curve_heuristics(intermed_result_struct.r1, ...
                                           intermed_result_struct.r2);
            
            % Get the figures created by pareto functions
            pareto_heuristics_fig = figs(6);
            pareto_front_with_design_images_fig = figs(4);
            pareto_front_LCOE_contours_fig = figs(3);
            pareto_constraint_activity_fig_damping = figs(1);
            pareto_constraint_activity_fig_reactive = figs(7);
            
            % Adjust figure position
            pareto_constraint_activity_fig_damping.Position  = [1 41 1536 844.8000];
            pareto_constraint_activity_fig_reactive.Position = [1 41 1536 844.8000];
            
            fig_array = [pareto_front_with_design_images_fig, ...
                        pareto_front_LCOE_contours_fig, ...
                        pareto_heuristics_fig, ...
                        pareto_constraint_activity_fig_damping, ...
                        pareto_constraint_activity_fig_reactive];
            
            tab_array_display = {};
            tab_array_latex = {};
            
            tab_firstrows = {};
            tab_colspecs = {};
            
            end_result_struct = struct();
end
