classdef LocationSensitivity < GenericAnalysis
    % LOCATIONSENSITIVITY Analysis class for location sensitivity tables
    %   Generates location sensitivity analysis table

    properties
        fig_names = {'location_1_power_matrix', 'location_2_power_matrix', ...
                     'location_3_power_matrix', 'location_4_power_matrix', ...
                     'location_probability_PDF'}
        tab_names = {'location_sensitivity'}
    end

    methods (Static)

        function intermed_result_struct = analysis_fcn(p, b)
            % Run location sensitivity analysis
            [flux, BW, storm, depths, locs, ...
             most_common_wave, X_opts, ...
             obj_opts, flags, figs] = location_sensitivity(p, b);

            % Store results for post-processing
            intermed_result_struct.p = p;
            intermed_result_struct.b = b;
            intermed_result_struct.flux = flux;
            intermed_result_struct.BW = BW;
            intermed_result_struct.storm = storm;
            intermed_result_struct.depths = depths;
            intermed_result_struct.locs = locs;
            intermed_result_struct.most_common_wave = most_common_wave;
            intermed_result_struct.X_opts = X_opts;
            intermed_result_struct.obj_opts = obj_opts;
            intermed_result_struct.flags = flags;
            intermed_result_struct.h_power_matrix = figs(1:end - 1);
            intermed_result_struct.h_probability_PDF = figs(end);
        end

        function [fig_array, ...
                 tab_array_display, ...
                 tab_array_latex, ...
                 end_result_struct] = post_process_fcn(intermed_result_struct)

            p = intermed_result_struct.p;
            b = intermed_result_struct.b;
            flux = intermed_result_struct.flux;
            BW = intermed_result_struct.BW;
            storm = intermed_result_struct.storm;
            depths = intermed_result_struct.depths;
            locs = intermed_result_struct.locs;
            most_common_wave = intermed_result_struct.most_common_wave;
            X_opts = intermed_result_struct.X_opts;
            obj_opts = intermed_result_struct.obj_opts;
            flags = intermed_result_struct.flags;

            [tab, pct_diff, ...
             location_flags, ...
             tablatex] = location_post_process(p, b, flux, BW, storm, depths, ...
                                               locs, most_common_wave, ...
                                               X_opts, obj_opts, flags);

            fig_array = [intermed_result_struct.h_power_matrix; ...
                         intermed_result_struct.h_probability_PDF];

            tab_array_display = {tab};
            tab_array_latex = {tablatex};

            end_result_struct.location_flags = location_flags;
            end_result_struct.pct_diff_lcoe_hawaii = pct_diff;
        end

    end
end
