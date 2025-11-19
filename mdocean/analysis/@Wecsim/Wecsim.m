classdef Wecsim < GenericAnalysis
    %WECSIM Analysis class for WEC-Sim validation figures and tables
    %   Generates WEC-Sim error histograms and validation tables

    properties
        fig_names = {'WECSim_error_histograms_singlebody', 'WECSim_error_histograms_multibody', 'wecsim_all_sea_states'};
        tab_names = {'WECSim_errors'};
    end

    methods (Static)
        intermed_result_struct = analysis_fcn(~,~)

        [fig_array,...
                 tab_array_display,...
                 tab_array_latex,...
                 end_result_struct] = post_process_fcn(intermed_result_struct)
    end
end
