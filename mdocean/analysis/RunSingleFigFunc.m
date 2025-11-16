classdef RunSingleFigFunc < GenericAnalysis
    %RUNSINGLEFIGFUNC Analysis class for single run figures
    %   Generates nominal geometry, power matrix, and drag convergence figures
    
    properties
        fig_names = {'nominal_geometry_viz',...
                     'nominal_power_pdf',...
                     'nominal_power_matrix',...
                     'drag_convergence',...
                     'slamming_model_comparison'};
        tab_names = {};
    end
    
    methods (Static)
        
        function intermed_result_struct = analysis_fcn(p,b)
            % Run single analysis with drag convergence plotting
            p.operational_dynamics_debug_plots_on = true;
            figs = run_single(p,b);

            % Store figure numbers for post-processing
            intermed_result_struct.created_figs = figs;
        end
        
        function [fig_array,...
                 tab_array_display,...
                 tab_array_latex,...
                 end_result_struct] = post_process_fcn(intermed_result_struct)

            fig_array = intermed_result_struct.created_figs;
            
            tab_array_display = {};
            tab_array_latex = {};
            
            end_result_struct.single_run_complete = true;
        end
        
    end
end
