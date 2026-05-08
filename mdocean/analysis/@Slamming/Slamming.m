classdef Slamming < GenericAnalysis
    %SLAMMING Analysis class for slamming figures
    %   Generates slamming amplitude figures

    properties
        fig_names = {'slamming_amplitude'};
        tab_names = {};
    end

    methods
        function obj = Slamming(varargin)
            obj = obj@GenericAnalysis(varargin{:});
            obj.extra_analysis_outputs = { ...
                'results/Slamming/intermed_figure_handle_1.fig' ...
            };
        end
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
