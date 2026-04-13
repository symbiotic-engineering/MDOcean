classdef FitOlaya < GenericAnalysis
    %FITOLAYA Analysis class for Olaya damping plate hydrodynamic coefficient fits
    %   Loads digitized data and creates derived signals in analysis_fcn;
    %   generates all fit and exploratory plots in post_process_fcn.

    properties
        fig_names = {};
        tab_names = {};
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
