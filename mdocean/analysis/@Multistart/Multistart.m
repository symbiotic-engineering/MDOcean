classdef Multistart < GenericAnalysis
    %MULTISTART Analysis class for multistart optimization figures and tables
    %   Generates multistart convergence tree and parallel coordinates figures with results table
    
    properties
        fig_names = {'multistart_convergence_tree',...
                     'multistart_parallel_coordinates',...
                     'multistart_bar_chart'};
        tab_names = {'multistart_results'};
    end
    
    methods (Static)
        
        intermed_result_struct = analysis_fcn(p,b)
        
        [fig_array,...
        tab_array_display,...
        tab_array_latex,...
        end_result_struct,...
        tab_firstrows,...
        tab_colspecs] = post_process_fcn(intermed_result_struct)
        
    end
end
