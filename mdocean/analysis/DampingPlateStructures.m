classdef DampingPlateStructures < GenericAnalysis
    %DAMPINGPLATESTRUCTURES Analysis class for damping plate structural figures
    %   Generates damping plate moment, deflection, and aspect ratio figures

    properties
        fig_names = {'damping_plate_aspect_ratio', 'damping_plate_deflection', 'damping_plate_moment'};
        tab_names = {};
    end

    methods (Static)

        function intermed_result_struct = analysis_fcn(~,~)
            % Run damping plate structural analysis
            addpath('../dev/structures/damping-plate');
            BoedoPrantilAnnularPlate()

            % Store figure numbers for post-processing
            n = gcf().Number;
            intermed_result_struct.figs = [figure(n), figure(n-1), figure(n-7)];
        end

        function [fig_array,...
                 tab_array_display,...
                 tab_array_latex,...
                 end_result_struct] = post_process_fcn(intermed_result_struct)

            fig_array = intermed_result_struct.figs;

            tab_array_display = {};
            tab_array_latex = {};

            end_result_struct.damping_plate_analysis_complete = true;
        end

    end
end
