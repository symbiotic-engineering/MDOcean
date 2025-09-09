classdef RunSingleFigFunc < GenericAnalysis
    %RUNSINGLEFIGFUNC Analysis class for single run figures
    %   Generates nominal geometry, power matrix, and drag convergence figures
    
    properties
        fig_names = {'nominal_geometry_viz', 'nominal_power_pdf', 'nominal_power_matrix', 'drag_convergence', 'slamming_model_comparison'};
        tab_names = {};
    end
    
    methods (Static)
        
        function intermed_result_struct = analysis_fcn(p,b)
            % Run single analysis with drag convergence plotting
            p.operational_drag_convergence_plot_on = true;
            run_single(p,b)

            % Store figure numbers for post-processing
            intermed_result_struct.final_figure_number = gcf().Number;
        end
        
        function [fig_array,...
                 tab_array_display,...
                 tab_array_latex,...
                 end_result_struct] = post_process_fcn(intermed_result_struct)
            
            n = intermed_result_struct.final_figure_number;
            
            fig_array = [figure(n), figure(n-1), figure(n-2), figure(n-3), figure(n-4)];
            
            tab_array_display = {};
            tab_array_latex = {};
            
            end_result_struct.single_run_complete = true;
        end
        
    end
end
