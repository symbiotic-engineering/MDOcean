classdef ParetoFigFunc < GenericAnalysis
    % PARETOFIGFUNC Analysis class for Pareto optimization figures
    %   Generates Pareto front figures and heuristics visualization

    properties
        fig_names = {'pareto_front_with_design_images', 'pareto_front_LCOE_contours', ...
                     'pareto_heuristics', 'pareto_constraint_activity_damping', ...
                     'pareto_constraint_activity_reactive'}
        tab_names = {}
    end

    methods (Static)

        function intermed_result_struct = analysis_fcn(p, b)
            % Run Pareto search
            [r1_damping, r2_reactive] = damping_vs_reactive(p, b);

            intermed_result_struct.r1 = r1_damping;
            intermed_result_struct.r2 = r2_reactive;
        end

        function [fig_array, ...
                 tab_array_display, ...
                 tab_array_latex, ...
                 end_result_struct] = post_process_fcn(intermed_result_struct)

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

            end_result_struct = struct();
        end

    end
end
