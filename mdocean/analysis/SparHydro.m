classdef SparHydro < GenericAnalysis
    %SPARHYDRO Analysis class for spar hydrodynamic figures
    %   Generates spar added mass figures
    
    properties
        fig_names = {'spar_added_mass'};
        tab_names = {};
    end
    
    methods (Static)
        
        function intermed_result_struct = analysis_fcn(~,~)
            % Run spar hydrodynamic analysis
            fig = spar_hydro_plot();
            
            % Store figure for post-processing
            intermed_result_struct.figure_handle = fig;
        end
        
        function [fig_array,...
                 tab_array_display,...
                 tab_array_latex,...
                 end_result_struct] = post_process_fcn(intermed_result_struct)
            
            fig_array = intermed_result_struct.figure_handle;
            
            tab_array_display = {};
            tab_array_latex = {};
            
            end_result_struct.spar_hydro_analysis_complete = true;
        end
        
    end
end
