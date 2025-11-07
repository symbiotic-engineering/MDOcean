classdef DescFcns < GenericAnalysis
    % DESCFCNS Analysis class for describing function figures
    %   Generates drag and saturation describing function figures

    properties
        fig_names = {'drag_desc_fcn', 'saturation_desc_fcn'}
        tab_names = {}
    end

    methods (Static)

        function intermed_result_struct = analysis_fcn(~, ~)
            % Run describing function demo and capture returned figure handles
            figs = sin_desc_fcn_demo();
            intermed_result_struct.created_figs = figs;
        end

        function [fig_array, ...
                 tab_array_display, ...
                 tab_array_latex, ...
                 end_result_struct] = post_process_fcn(intermed_result_struct)

            % Return the figures created by the demo; GenericAnalysis will
            % validate and map them to fig_names.
            fig_array = intermed_result_struct.created_figs;

            tab_array_display = {};
            tab_array_latex = {};

            end_result_struct.desc_fcn_analysis_complete = true;
        end

    end
end
