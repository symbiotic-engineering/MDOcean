classdef ExampleAnalysis < GenericAnalysis
    % EXAMPLEANALYSIS Summary of this class goes here
    %   Detailed explanation goes here

    properties
        fig_names = {'My_fig_1', 'My_fig_2'}
        tab_names = {'My_tab'}
    end

    methods (Static)

        function intermed_result_struct = analysis_fcn()
            disp('doing time consuming analysis');
            pause(2);
            intermed_result_struct.x = [1 2 3];
            intermed_result_struct.y = [4 5 6];
        end

        function [fig_array, ...
         tab_array_display, ...
         tab_array_latex, ...
         end_result_struct] = post_process_fcn(intermed_result_struct)
            f1 = figure;
            plot(intermed_result_struct.x, intermed_result_struct.y);
            f2 = figure;
            plot(intermed_result_struct.x, -intermed_result_struct.y);
            fig_array = [f1, f2];

            tab_array_display = {table(intermed_result_struct.x, intermed_result_struct.y)};
            tab_array_latex = tab_array_display;

            end_result_struct.z = intermed_result_struct.x + intermed_result_struct.y;
        end

    end
end
