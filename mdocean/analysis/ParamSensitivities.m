classdef ParamSensitivities < GenericAnalysis
    %PARAMSENSITIVITIES Analysis class for parameter sensitivity figures
    %   Generates parameter sensitivity analysis figures including tornado plots and grids
    
    properties
        fig_names = {'re_optim_constraint_grid', 'post_optim_constraint_grid', ...
                    're_optim_design_grid', 'post_optim_design_grid', ...
                    're_optim_objective_grid', 'post_optim_re_optim_objective_grid', ...
                    're_optim_design_tornado_J2', 're_optim_design_tornado_J1', ...
                    're_optim_objective_tornado', 'nonlinear_design_J2', ...
                    'nonlinear_design_J1', 'nonlinear_objectives'};
        tab_names = {};
    end
    
    methods (Static)
        
        function intermed_result_struct = analysis_fcn(p,b)
            % Run parameter sweep analysis
            [fig_arr, runtime_post_optim, runtime_re_optim] = param_sweep(b.filename_uuid);
            
            % Store results for post-processing
            intermed_result_struct.fig_arr = fig_arr;
            intermed_result_struct.runtime_post_optim = runtime_post_optim;
            intermed_result_struct.runtime_re_optim = runtime_re_optim;
            intermed_result_struct.num_DVs = length(b.var_names) - 1;
            intermed_result_struct.var_names = b.var_names;
        end
        
        function [fig_array,...
                 tab_array_display,...
                 tab_array_latex,...
                 end_result_struct] = post_process_fcn(intermed_result_struct)
            
            fig_arr = intermed_result_struct.fig_arr;
            num_DVs = intermed_result_struct.num_DVs;
            var_names = intermed_result_struct.var_names;
            
            % Extract figures in the expected order
            fig_array = [];
            
            % Main figures
            fig_array = [fig_array, fig_arr(end-5)]; % post_optim_re_optim_objective_grid
            fig_array = [fig_array, fig_arr(end-4)]; % re_optim_objective_grid
            fig_array = [fig_array, fig_arr(end-3)]; % post_optim_design_grid
            fig_array = [fig_array, fig_arr(end-2)]; % re_optim_design_grid
            fig_array = [fig_array, fig_arr(end-1)]; % post_optim_constraint_grid
            fig_array = [fig_array, fig_arr(end)];   % re_optim_constraint_grid
            
            % Tornado plots
            fig_array = [fig_array, fig_arr(end-7)]; % re_optim_design_tornado_J1
            fig_array = [fig_array, fig_arr(end-6)]; % re_optim_design_tornado_J2
            
            % Design variable specific tornado plots
            for i = 1:num_DVs
                fig_array = [fig_array, fig_arr(end-(7+i))];
            end
            
            % Additional plots
            fig_array = [fig_array, fig_arr(end-num_DVs-8)];  % re_optim_objective_tornado
            fig_array = [fig_array, fig_arr(end-num_DVs-9)];  % nonlinear_design_J2
            fig_array = [fig_array, fig_arr(end-num_DVs-10)]; % nonlinear_design_J1
            fig_array = [fig_array, fig_arr(end-num_DVs-11)]; % nonlinear_objectives
            
            tab_array_display = {};
            tab_array_latex = {};
            
            end_result_struct.sensitivity_analysis_complete = true;
            end_result_struct.runtime_post_optim = intermed_result_struct.runtime_post_optim;
            end_result_struct.runtime_re_optim = intermed_result_struct.runtime_re_optim;
        end
        
    end
end
