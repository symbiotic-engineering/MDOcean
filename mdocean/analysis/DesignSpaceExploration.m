classdef DesignSpaceExploration < GenericAnalysis
    %DESIGNSPACEEXPLORATION Analysis class for design space exploration figures
    %   Generates design space exploration experimental figures

    properties
        fig_names = {'experiments_pareto','experiments_ratios'};
        tab_names = {'experiments_results'};
    end

    methods (Static)

        function intermed_result_struct = analysis_fcn(p,b)
            % Run design space exploration experiments
            [X_ins,ratios,LCOE,cost,power,failed] = experiments(p,b);

            cd('..')
            pareto = ParetoFigFunc(p,b);
            pareto = pareto.run_analysis_from_load_if_possible();
            pareto_results_struct = pareto.intermed_result_struct.r1;
            cd('mdocean')

            % Store figure for post-processing
            intermed_result_struct.b = b;
            intermed_result_struct.X_ins = X_ins;
            intermed_result_struct.ratios = ratios;
            intermed_result_struct.LCOE = LCOE;
            intermed_result_struct.cost = cost;
            intermed_result_struct.power = power;
            intermed_result_struct.failed = failed;
            intermed_result_struct.pareto_results_struct = pareto_results_struct;
        end

        function [fig_array,...
                 tab_array_display,...
                 tab_array_latex,...
                 end_result_struct] = post_process_fcn(intermed_result_struct)

            b = intermed_result_struct.b;
            X_ins = intermed_result_struct.X_ins;
            ratios = intermed_result_struct.ratios;
            LCOE = intermed_result_struct.LCOE;
            cost = intermed_result_struct.cost;
            power = intermed_result_struct.power;
            failed = intermed_result_struct.failed;
            pareto_results_struct = intermed_result_struct.pareto_results_struct;

            [fig_array,results_tab] = experiments_plot(b,X_ins,ratios,LCOE,...
                                                        cost,power,failed,...
                                                        pareto_results_struct);

            tab_array_display = {results_tab};
            tab_array_latex = {};

            end_result_struct.design_space_exploration_complete = true;
        end

    end
end
