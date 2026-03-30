classdef SparHydro < GenericAnalysis
    %SPARHYDRO Analysis class for spar hydrodynamic figures
    %   Generates spar added mass figures

    properties
        fig_names = {'spar_added_mass'};
        tab_names = {};
    end

    methods
        function obj = SparHydro(varargin)
            obj = obj@GenericAnalysis(varargin{:});
            obj.extra_outputs = { ...
                'results/SparHydro/intermed_figure_handle_1.fig' ...
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
