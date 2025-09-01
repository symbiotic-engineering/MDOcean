classdef (Abstract) GenericAnalysis
    %GENERICANALYSIS Summary of this class goes here
    %   Detailed explanation goes here

    properties
        p = parameters();
        b = var_bounds();
        fig_array
        tab_array_display
        tab_array_latex
        intermed_result_struct % all results needed to debug and make plots
        end_result_struct % only summary end results for latex output or tests
        output_folder
        analysis_dependencies
        postpro_dependencies
        analysis_outputs
        postpro_outputs
    end

    properties (Abstract)
        fig_names
        tab_names
    end

    methods (Abstract, Static)

        intermed_result_struct = analysis_fcn();

        [fig_array,...
         tab_array_display,...
         tab_array_latex,...
         end_result_struct] = post_process_fcn(intermed_result_struct);

    end

    methods
        function obj = GenericAnalysis()
            obj.output_folder = ['results' filesep class(obj)];
            mkdir(obj.output_folder)

            obj.analysis_dependencies = matlab.codetools.requiredFilesAndProducts([class(obj) '.analysis_fcn()']);
            obj.analysis_outputs = {[obj.output_folder, filesep, 'intermed.mat']};
            obj.postpro_dependencies = [matlab.codetools.requiredFilesAndProducts([class(obj) '.post_process_fcn()']),...
                                        obj.analysis_outputs];

            obj.postpro_outputs  = strcat(obj.output_folder,...
                filesep, [strcat(obj.fig_names,'.pdf') ...
                          strcat(obj.tab_names,'.tex') ...
                          'end.mat']);
        end

        function obj = run_analysis(obj)
            obj.intermed_result_struct = obj.analysis_fcn();
            obj.save_intermed_results();
        end

        function obj = run_post_process(obj)
            [obj.fig_array,...
             obj.tab_array_display,...
             obj.tab_array_latex,...
             obj.end_result_struct] = obj.post_process_fcn(obj.intermed_result_struct);

            obj.save_figs();
            obj.save_tabs();
            obj.save_end_results();
        end

        function obj = load_intermed_results(obj)
            obj.intermed_result_struct = load([obj.output_folder filesep 'intermed']);
        end

        function save_intermed_results(obj)
            s = obj.intermed_result_struct;
            save([obj.output_folder filesep 'intermed'], '-struct', 's')
        end

        function save_end_results(obj)
            s =  obj.end_result_struct;
            save([obj.output_folder filesep 'end'],'-struct', 's')
        end

        function save_figs(obj)
            for i=1:length(obj.fig_names)
                fig = obj.fig_array(i);
                fname = [obj.output_folder filesep obj.fig_names{i}];
                save_pdf(fig, fname)
            end
        end
        
        function save_tabs(obj)
            for i=1:length(obj.tab_names)
                tab = obj.tab_array_latex{i};
                fname = [obj.output_folder filesep obj.tab_names{i}];
                table2latex(tab,fname)
            end
        end

        function run_all_from_analysis(obj)
            obj = run_analysis(obj);
            obj.run_post_process()
        end

        function run_all_from_load(obj)
            obj.load_intermed_results();
            obj.run_post_process();
        end

        function write_calkit(obj)
            cell2filelist = @(c) char(join(strcat("    - ",c),newline));

            analysis_stage = ['stage: analysis-' class(obj) newline ...
                              '  kind: matlab-command' newline ...
                              '  environment: matlab' newline ...
                              '  command: obj=' class(obj) '; obj.run_analysis();' newline ...
                              '  inputs: ' newline ...
                              cell2filelist(obj.analysis_dependencies.') newline ...
                              '  outputs: ' newline ...
                              cell2filelist(obj.analysis_outputs.') ];
            postpro_stage = ['stage: postpro-' class(obj) newline ...
                              '  kind: matlab-command' newline ...
                              '  environment: matlab' newline ...
                              '  command: obj=' class(obj) '; obj.run_all_from_load();' newline ...
                              '  inputs: ' newline ...
                              cell2filelist(obj.postpro_dependencies.') newline ...
                              '  outputs: ' newline ...
                              cell2filelist(obj.postpro_outputs.') ];
            stages = [analysis_stage newline postpro_stage];

            disp(stages)
        end
    end
end