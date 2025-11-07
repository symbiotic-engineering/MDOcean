classdef JpdMultiplyFigFunc < GenericAnalysis
    % JPDMULTIPLYFIGFUNC Analysis class for JPD multiplication figures
    %   Generates JPD multiplication power matrix figure

    properties
        fig_names = {'JPD_multiplication'}
        tab_names = {}
    end

    methods (Static)

        function intermed_result_struct = analysis_fcn(p, b)
            % Generate JPD multiplication figure
            X = [b.X_noms; 1];

            h = plot_power_matrix(X, p, b, b.filename_uuid);
            intermed_result_struct.created_figs = h;
        end

        function [fig_array, ...
                 tab_array_display, ...
                 tab_array_latex, ...
                 end_result_struct] = post_process_fcn(intermed_result_struct)

            fig_array = intermed_result_struct.created_figs;

            tab_array_display = {};
            tab_array_latex = {};

            end_result_struct.jpd_multiplication_complete = true;
        end

    end
end
