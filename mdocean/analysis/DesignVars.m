classdef DesignVars < GenericAnalysis
    %DESIGNVARS Analysis class for design variables table
    %   Generates design variables bounds table

    properties
        fig_names = {};
        tab_names = {'design_vars'};
    end

    methods (Static)

        function intermed_result_struct = analysis_fcn(~,b)
            % Get design variables information
            % Store results for post-processing
            intermed_result_struct.X_mins = b.X_mins;
            intermed_result_struct.X_noms = b.X_noms;
            intermed_result_struct.X_maxs = b.X_maxs;
            intermed_result_struct.var_names = b.var_names;
        end

        function [fig_array,...
                 tab_array_display,...
                 tab_array_latex,...
                 end_result_struct] = post_process_fcn(intermed_result_struct)

            tab = array2table([intermed_result_struct.X_mins intermed_result_struct.X_noms intermed_result_struct.X_maxs], ...
                    'VariableNames',{'Mins','Noms','Maxs'}, 'RowNames', intermed_result_struct.var_names(1:end-1));

            fig_array = [];

            tab_array_display = {tab};
            tab_array_latex = {tab};

            end_result_struct.design_vars_table = tab;
        end

    end
end
