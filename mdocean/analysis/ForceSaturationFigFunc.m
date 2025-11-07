classdef ForceSaturationFigFunc < GenericAnalysis
    %FORCESATURATIONFIGFUNC Analysis class for force saturation figures
    %   Generates power-force sensitivity and runtime sensitivity figures

    properties
        fig_names = {'power_force_sensitivity', 'runtime_sensitivity'};
        tab_names = {};
    end

    methods (Static)

        function intermed_result_struct = analysis_fcn(p,b)
            % Run force saturation analysis
            [fig1, fig2] = force_sat_results(p, b);

            % Store figures for post-processing
            intermed_result_struct.fig1 = fig1;
            intermed_result_struct.fig2 = fig2;
        end

        function [fig_array,...
                 tab_array_display,...
                 tab_array_latex,...
                 end_result_struct] = post_process_fcn(intermed_result_struct)

            fig_array = [intermed_result_struct.fig1, intermed_result_struct.fig2];

            tab_array_display = {};
            tab_array_latex = {};

            end_result_struct.force_saturation_analysis_complete = true;
        end

    end
end
