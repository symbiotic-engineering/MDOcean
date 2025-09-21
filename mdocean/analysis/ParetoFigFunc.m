classdef ParetoFigFunc < GenericAnalysis
    %PARETOFIGFUNC Analysis class for Pareto optimization figures
    %   Generates Pareto front figures and heuristics visualization
    
    properties
        fig_names = {'pareto_front_with_design_images', 'pareto_front_LCOE_contours', ...
                    'pareto_heuristics', 'pareto_constraint_activity'};
        tab_names = {};
    end
    
    methods (Static)
        
        function intermed_result_struct = analysis_fcn(~,b)
            % Run Pareto search and heuristics
            pareto_search(b.filename_uuid)
            pareto_curve_heuristics()
            
            % Store figure numbers for post-processing
            intermed_result_struct.heuristics_fig_number = gcf().Number;
        end
        
        function [fig_array,...
                 tab_array_display,...
                 tab_array_latex,...
                 end_result_struct] = post_process_fcn(intermed_result_struct)
            
            n = intermed_result_struct.heuristics_fig_number;
            
            % Get the figures created by pareto functions
            pareto_heuristics_fig = figure(n);
            pareto_front_with_design_images_fig = figure(n - 2);
            pareto_front_LCOE_contours_fig = figure(n - 3);
            pareto_constraint_activity_fig = figure(n - 5);
            
            % Adjust figure position
            pareto_constraint_activity_fig.Position = [1 41 1536 844.8000];
            
            fig_array = [pareto_front_with_design_images_fig, pareto_front_LCOE_contours_fig, ...
                        pareto_heuristics_fig, pareto_constraint_activity_fig];
            
            tab_array_display = {};
            tab_array_latex = {};
            
            end_result_struct.pareto_points_generated = true;
        end
        
    end
end
