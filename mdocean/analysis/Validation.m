classdef Validation < GenericAnalysis
    %VALIDATION Analysis class for validation figures and tables
    %   Generates cost vs N WEC validation figure and validation table
    
    properties
        fig_names = {'cost_vs_N_WEC'};
        tab_names = {'validation'};
    end
    
    methods (Static)
        
        function intermed_result_struct = analysis_fcn()
            % Run validation analysis
            [~, ~, ~, ~, tab1a, fig_cost_vs_N_WEC] = validate_nominal_RM3('report');
            [~,~,~,~,tab1b] = validate_nominal_RM3('wecsim');
            
            % Store results for post-processing
            intermed_result_struct.tab1a = tab1a;
            intermed_result_struct.tab1b = tab1b;
            intermed_result_struct.fig_cost_vs_N_WEC = fig_cost_vs_N_WEC;
        end
        
        function [fig_array,...
                 tab_array_display,...
                 tab_array_latex,...
                 end_result_struct] = post_process_fcn(intermed_result_struct)
            
            tab1a = intermed_result_struct.tab1a;
            tab1b = intermed_result_struct.tab1b;
            
            % merge table 1a and 1b while preserving row order
            sharedCols = intersect(tab1a.Properties.VariableNames, tab1b.Properties.VariableNames);
            tab1a.RowNum = (1:length(tab1a.Properties.RowNames))';
            tab1b.RowNum = length(tab1a.Properties.RowNames) + (1:length(tab1b.Properties.RowNames))';
            
            tab1 = outerjoin(tab1a, tab1b, 'Keys', [{'RowNum'},sharedCols], 'MergeKeys', true);
            tab1 = removevars(tab1,'RowNum');
            tab1.Properties.RowNames = [tab1a.Properties.RowNames; tab1b.Properties.RowNames];

            vector_cols = {'capex','capex_design','J_capex_design','capex_struct','capex_PTO','opex','LCOE'};
            idx_remove = ismember(tab1.Properties.VariableNames,vector_cols); 
            tab1latex = rows2vars(tab1,'VariableNamingRule','preserve','DataVariables',~idx_remove);
            tab1latex.OriginalVariableNames = remove_underscores(modify_suffix(tab1latex.OriginalVariableNames));
            new_names = {'Variable','MDOcean','Actual','Error','MDOcean ','Actual ','Error '};
            tab1latex = renamevars(tab1latex, tab1latex.Properties.VariableNames, new_names);
            
            fig_array = intermed_result_struct.fig_cost_vs_N_WEC;
            
            tab_array_display = {tab1};
            tab_array_latex = {tab1latex};
            
            end_result_struct.validation_complete = true;
            end_result_struct.validation_table = tab1;
        end
        
    end
end
