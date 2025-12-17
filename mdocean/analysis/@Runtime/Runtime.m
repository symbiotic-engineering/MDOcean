classdef Runtime < GenericAnalysis
    %RUNTIME Analysis class for runtime comparison figures
    %   Generates dynamics, hydro, and simulation runtime figures

    properties
        fig_names = {'dynamics_runtime', 'hydro_runtime', 'sim_runtime'};
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
