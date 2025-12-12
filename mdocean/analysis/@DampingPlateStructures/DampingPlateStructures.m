classdef DampingPlateStructures < GenericAnalysis
    %DAMPINGPLATESTRUCTURES Analysis class for damping plate structural figures
    %   Generates damping plate moment, deflection, and aspect ratio figures

    properties
        fig_names = {'damping_plate_aspect_ratio', 'damping_plate_deflection', 'damping_plate_moment'};
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
