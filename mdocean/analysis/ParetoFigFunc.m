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
            [x,fval,residuals,tol,p] = pareto_search(b.filename_uuid);
            
            % Store results
            intermed_result_struct.x = x;
            intermed_result_struct.fval = fval;
            intermed_result_struct.residuals = residuals;
            intermed_result_struct.tol = tol;
            intermed_result_struct.p = p;
        end
        
        function [fig_array,...
                 tab_array_display,...
                 tab_array_latex,...
                 end_result_struct] = post_process_fcn(intermed_result_struct)
            
            figs = pareto_curve_heuristics(intermed_result_struct);
            
            % Get the figures created by pareto functions
            pareto_heuristics_fig = figs(6);
            pareto_front_with_design_images_fig = figs(4);
            pareto_front_LCOE_contours_fig = figs(3);
            pareto_constraint_activity_fig = figs(1);
            
            % Adjust figure position
            pareto_constraint_activity_fig.Position = [1 41 1536 844.8000];
            
            fig_array = [pareto_front_with_design_images_fig, ...
                        pareto_front_LCOE_contours_fig, ...
                        pareto_heuristics_fig, ...
                        pareto_constraint_activity_fig];
            
            tab_array_display = {};
            tab_array_latex = {};
            
            end_result_struct = struct();
        end
        
    end
end
