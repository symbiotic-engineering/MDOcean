classdef DesignSpaceExploration < GenericAnalysis
    %DESIGNSPACEEXPLORATION Analysis class for design space exploration figures
    %   Generates design space exploration experimental figures
    
    properties
        fig_names = {'experiments'};
        tab_names = {};
    end
    
    methods (Static)
        
        function intermed_result_struct = analysis_fcn(p,b)
            % Run design space exploration experiments and capture the
            % figures returned by the experiments helper. experiments()
            % already returns handles when available.
            try
                figs = experiments(p,b);
            catch err
                warning('DesignSpaceExploration:experimentsCall','Failed to run experiments(): %s', err.message);
                figs = [];
            end

            if iscell(figs)
                figs = [figs{:}];
            end

            intermed_result_struct.created_figs = figs;
        end
        
        function [fig_array,...
                 tab_array_display,...
                 tab_array_latex,...
                 end_result_struct] = post_process_fcn(intermed_result_struct)
            
            fig_array = intermed_result_struct.created_figs;
            
            tab_array_display = {};
            tab_array_latex = {};
            
            end_result_struct.design_space_exploration_complete = true;
        end
        
    end
end
