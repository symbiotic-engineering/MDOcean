classdef Runtime < GenericAnalysis
    %RUNTIME Analysis class for runtime comparison figures
    %   Generates dynamics, hydro, and simulation runtime figures
    
    properties
        fig_names = {'dynamics_runtime', 'hydro_runtime', 'sim_runtime'};
        tab_names = {};
    end
    
    methods
        
        function intermed_result_struct = analysis_fcn(obj)
            % Run runtime comparison analysis
            [f1, f2, f3] = module_runtime_compare(obj.p,obj.b);
            
            % Store results for post-processing
            intermed_result_struct.dynamics_runtime_fig = f1;
            intermed_result_struct.hydro_runtime_fig = f2;
            intermed_result_struct.sim_runtime_fig = f3;
        end
        
        function [fig_array,...
                 tab_array_display,...
                 tab_array_latex,...
                 end_result_struct] = post_process_fcn(intermed_result_struct)
            
            fig_array = [intermed_result_struct.dynamics_runtime_fig, ...
                        intermed_result_struct.hydro_runtime_fig, ...
                        intermed_result_struct.sim_runtime_fig];
            
            tab_array_display = {};
            tab_array_latex = {};
            
            end_result_struct.runtime_analysis_complete = true;
        end
        
    end
end
