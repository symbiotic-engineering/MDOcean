classdef Comparison < GenericAnalysis
    %COMPARISON Analysis class for comparison figures and tables
    %   Generates overlaid geometry, hydro coeffs, probability CDF figures and optimal design tables
    
    properties
        fig_names = {'overlaid_geometry', 'overlaid_hydro_coeffs', 'probability_CDF'};
        tab_names = {'optimal_design_vars', 'optimal_outputs'};
    end
    
    methods (Static)
        
        function intermed_result_struct = analysis_fcn(p,b)
            % Run comparison analysis
            [optimal_design_vars, optimal_outputs] = compare(p,b);
            
            % Store results for post-processing
            intermed_result_struct.optimal_design_vars = optimal_design_vars;
            intermed_result_struct.optimal_outputs = optimal_outputs;
            intermed_result_struct.final_figure_number = gcf().Number;
        end
        
        function [fig_array,...
                 tab_array_display,...
                 tab_array_latex,...
                 end_result_struct] = post_process_fcn(intermed_result_struct)
            
            n = intermed_result_struct.final_figure_number;
            
            % Get figures in the expected order
            fig_array = [figure(n-3), figure(n), figure(n-2)];
            
            tab_array_display = {intermed_result_struct.optimal_design_vars, ...
                               intermed_result_struct.optimal_outputs};
            tab_array_latex = tab_array_display;
            
            end_result_struct.comparison_complete = true;
            end_result_struct.optimal_design_vars = intermed_result_struct.optimal_design_vars;
            end_result_struct.optimal_outputs = intermed_result_struct.optimal_outputs;
        end
        
    end
end
