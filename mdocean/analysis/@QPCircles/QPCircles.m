classdef QPCircles < GenericAnalysis
    %QPCIRCLES Analysis class for QP circle intersection figure
    %   Generates the QP circles figure from circle_intersect_script

    properties
        fig_names = {'qp_circles'};
        tab_names = {};
    end
    methods
        function obj = QPCircles(varargin)
            obj = obj@GenericAnalysis(varargin{:});
            obj.extra_analysis_outputs = { ...
                'results/QPCircles/intermed_figure_handle_1.fig'};
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
