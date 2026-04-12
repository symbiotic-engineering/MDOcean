classdef PtoSweep < GenericAnalysis
    %PTOSWEEP Analysis class for PTO parameter sweep figures
    %   Performs a 2-D sweep of F_max (maximum powertrain force) and P_max
    %   (maximum powertrain power) for the nominal RM3 geometry.  For each
    %   grid point the full simulation is evaluated and the resulting
    %   average power, design cost, and LCOE are stored.  The
    %   post-processing function produces a single figure with three
    %   contourf panels (power, cost, LCOE) overlaid with a hatched
    %   infeasibility region and markers for the maximum-power and
    %   minimum-LCOE operating points.

    properties
        fig_names = {'pto_sweep'};
        tab_names = {};
    end

    methods (Static)
        intermed_result_struct = analysis_fcn(p, b)

        [fig_array, ...
         tab_array_display, ...
         tab_array_latex, ...
         end_result_struct, ...
         tab_firstrows, ...
         tab_colspecs] = post_process_fcn(intermed_result_struct)
    end
end
