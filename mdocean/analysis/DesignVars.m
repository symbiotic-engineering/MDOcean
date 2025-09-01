classdef DesignVars < GenericAnalysis
    %DESIGNVARS Analysis class for design variables table
    %   Generates design variables bounds table
    
    properties
        fig_names = {};
        tab_names = {'design_vars'};
    end
    
    methods
        
        function intermed_result_struct = analysis_fcn(obj)
            % Get design variables information
            % Store results for post-processing
            intermed_result_struct.X_mins = obj.b.X_mins;
            intermed_result_struct.X_noms = obj.b.X_noms;
            intermed_result_struct.X_maxs = obj.b.X_maxs;
            intermed_result_struct.var_names = obj.b.var_names;
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
