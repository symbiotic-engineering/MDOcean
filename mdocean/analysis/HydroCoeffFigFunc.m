classdef HydroCoeffFigFunc < GenericAnalysis
    % HYDROCOEFFIGFUNC Analysis class for hydrodynamic coefficient figures
    %   Generates hydrodynamic coefficient error figures

    properties
        fig_names = {'hydro_coeff_err'}
        tab_names = {}
    end

    methods (Static)

        function intermed_result_struct = analysis_fcn(~, ~)
            % Run hydrodynamic coefficient error analysis
            [~, ~, fig] = hydro_coeff_err();

            % Store figure for post-processing
            intermed_result_struct.created_figs = fig;
        end

        function [fig_array, ...
                 tab_array_display, ...
                 tab_array_latex, ...
                 end_result_struct] = post_process_fcn(intermed_result_struct)

            % Return created figures (validation happens in GenericAnalysis)
            fig_array = intermed_result_struct.created_figs;

            tab_array_display = {};
            tab_array_latex = {};

            end_result_struct.hydro_coeff_analysis_complete = true;
        end

    end
end
