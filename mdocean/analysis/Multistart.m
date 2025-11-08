classdef Multistart < GenericAnalysis
    %MULTISTART Analysis class for multistart optimization figures and tables
    %   Generates multistart convergence tree and parallel coordinates figures with results table
    
    properties
        fig_names = {'multistart_convergence_tree', 'multistart_parallel_coordinates','multistart_bar_chart'};
        tab_names = {'multistart_results'};
    end
    
    methods (Static)
        
        function intermed_result_struct = analysis_fcn(p,b)
            % Run multistart optimization analysis
            [X_opt,objs,flags,x0s] = gradient_mult_x0(p,b);
            
            % Store results for post-processing
            intermed_result_struct.p = p;
            intermed_result_struct.b = b;
            intermed_result_struct.X_opt = X_opt;
            intermed_result_struct.objs = objs;
            intermed_result_struct.flags = flags;
            intermed_result_struct.x0s = x0s;
        end
        
        function [fig_array,...
                 tab_array_display,...
                 tab_array_latex,...
                 end_result_struct] = post_process_fcn(intermed_result_struct)
            
            p = intermed_result_struct.p;
            b = intermed_result_struct.b;
            X_opt =intermed_result_struct.X_opt;
            objs = intermed_result_struct.objs;
            flags = intermed_result_struct.flags;
            x0s = intermed_result_struct.x0s;

            [treeFig, parallelFig, barFig, results] = multistart_postpro(p,b,X_opt,objs,flags,x0s);

            fig_array = [treeFig, parallelFig, barFig];
            
            tab_array_display = {results};
            tab_array_latex = {results};

            end_result_struct = struct();
        end
        
    end
end
