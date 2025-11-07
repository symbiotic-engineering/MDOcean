classdef Wecsim < GenericAnalysis
    %WECSIM Analysis class for WEC-Sim validation figures and tables
    %   Generates WEC-Sim error histograms and validation tables

    properties
        fig_names = {'WECSim_error_histograms_singlebody', 'WECSim_error_histograms_multibody', 'wecsim_all_sea_states'};
        tab_names = {'WECSim_errors'};
    end

    methods (Static)

        function intermed_result_struct = analysis_fcn(~,~)
            % Run WEC-Sim dynamics validation
            [~, ~, ~, tab, fig_singlebody, fig_multibody] = validate_dynamics();

            % Store results for post-processing
            intermed_result_struct.validation_table = tab;
            intermed_result_struct.fig_singlebody = fig_singlebody;
            intermed_result_struct.fig_multibody = fig_multibody;
            intermed_result_struct.all_sea_states_fig = gcf();
        end

        function [fig_array,...
                 tab_array_display,...
                 tab_array_latex,...
                 end_result_struct] = post_process_fcn(intermed_result_struct)

            fig_array = [intermed_result_struct.fig_singlebody, ...
                        intermed_result_struct.fig_multibody, ...
                        intermed_result_struct.all_sea_states_fig];

            tab_array_display = {intermed_result_struct.validation_table};
            tab_array_latex = {intermed_result_struct.validation_table};

            end_result_struct.wecsim_validation_complete = true;
            end_result_struct.validation_errors = intermed_result_struct.validation_table;
        end

    end
end
