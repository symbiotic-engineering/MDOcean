classdef ForceSaturationFigFunc < GenericAnalysis
    %FORCESATURATIONFIGFUNC Analysis class for force saturation figures
    %   Generates power-force sensitivity and runtime sensitivity figures

    properties
        fig_names = {'power_force_sensitivity', 'runtime_sensitivity'};
        tab_names = {};
    end

    methods
        function obj = ForceSaturationFigFunc(varargin)
            obj = obj@GenericAnalysis(varargin{:});
            obj.extra_outputs = { ...
                'results/ForceSaturationFigFunc/intermed_fig2_1.fig', ...
                'results/ForceSaturationFigFunc/intermed_fig1_1.fig' ...
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
