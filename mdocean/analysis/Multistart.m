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
            [treeFig, parallelFig, barFig, tab] = gradient_mult_x0(p,b);

            % Store results for post-processing
            intermed_result_struct.treeFig = treeFig;
            intermed_result_struct.parallelFig = parallelFig;
            intermed_result_struct.barFig = barFig;
            intermed_result_struct.multistart_table = tab;
        end

        function [fig_array,...
                 tab_array_display,...
                 tab_array_latex,...
                 end_result_struct] = post_process_fcn(intermed_result_struct)

            fig_array = [intermed_result_struct.treeFig, intermed_result_struct.parallelFig, intermed_result_struct.barFig];

            tab_array_display = {intermed_result_struct.multistart_table};
            tab_array_latex = {intermed_result_struct.multistart_table};

            end_result_struct.multistart_complete = true;
            end_result_struct.multistart_results = intermed_result_struct.multistart_table;
        end

    end
end
