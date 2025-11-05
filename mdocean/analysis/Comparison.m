classdef Comparison < GenericAnalysis
    %COMPARISON Analysis class for comparison figures and tables
    %   Generates overlaid geometry, hydro coeffs, probability CDF figures and optimal design tables
    
    properties
        fig_names = {'overlaid_geometry', 'overlaid_hydro_coeffs', 'probability_CDF', 'comparison_power_matrix'};
        tab_names = {'optimal_design_vars', 'optimal_outputs'};
    end
    
    methods (Static)
        
        function intermed_result_struct = analysis_fcn(p,b)
            % Run comparison analysis
            [Xs,vals] = compare_run(p,b);

            % Store results and returned figure handles for post-processing
            intermed_result_struct.Xs = Xs;
            intermed_result_struct.vals = vals;
        end
        
        function [fig_array,...
                 tab_array_display,...
                 tab_array_latex,...
                 end_result_struct] = post_process_fcn(intermed_result_struct)
            
            Xs = intermed_result_struct.Xs;
            vals = intermed_result_struct.vals;
            [DV_table, out_table, fig_array] = compare(p,b,Xs,vals);
            
            tab_array_display = {DV_table, ...
                                 out_table};
            tab_array_latex = tab_array_display;
            
            end_result_struct.comparison_complete = true;
        end
        
    end
end
