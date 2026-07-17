classdef (Abstract) GroupedStageClass < handle
    %GROUPEDSTAGECLASS Wrapper exposing analysis + postprocessing as one stage.

    properties
        p
        b
        fig_array
        tab_array_display
        tab_array_latex
        tab_firstrows
        tab_colspecs
        intermed_result_struct
        end_result_struct
        output_folder
        analysis_outputs
        postpro_outputs
        extra_analysis_inputs
        extra_analysis_outputs
        extra_postpro_inputs
        extra_postpro_outputs

        analysis_stage
        postpro_stage
    end

    properties (Dependent)
        class_dependencies
        analysis_dependencies
        postpro_dependencies
    end

    properties (Abstract)
        fig_names
        tab_names
    end

    methods (Abstract, Static)
        intermed_result_struct = analysis_fcn(p,b);

        [fig_array,...
         tab_array_display,...
         tab_array_latex,...
         end_result_struct,...
         tab_firstrows,...
         tab_colspecs] = post_process_fcn(intermed_result_struct);
    end

    methods
        function obj = GroupedStageClass(p,b)
            if nargin < 1
                p = [];
            end
            if nargin < 2
                b = [];
            end
            class_name = class(obj);
            analysis_handle = str2func([class_name '.analysis_fcn']);
            postpro_handle = str2func([class_name '.post_process_fcn']);
            obj.analysis_stage = AnalysisClass(class_name, analysis_handle, p, b);
            if isempty(obj.tab_names)
                obj.postpro_stage = PostproClass(class_name, postpro_handle, obj.fig_names, obj.tab_names, obj.analysis_stage.p, obj.analysis_stage.b);
            else
                obj.postpro_stage = PostproTableClass(class_name, postpro_handle, obj.fig_names, obj.tab_names, obj.analysis_stage.p, obj.analysis_stage.b);
            end
            obj.copy_from_stages();
        end

        function val = get.class_dependencies(obj)
            val = unique([obj.analysis_stage.class_dependencies, obj.postpro_stage.class_dependencies]);
        end

        function val = get.analysis_dependencies(obj)
            obj.apply_to_stages();
            val = obj.analysis_stage.analysis_dependencies;
        end

        function val = get.postpro_dependencies(obj)
            obj.apply_to_stages();
            val = obj.postpro_stage.postpro_dependencies;
        end

        function obj = apply_to_stages(obj)
            obj.analysis_stage.p = obj.p;
            obj.analysis_stage.b = obj.b;
            obj.analysis_stage.output_folder = obj.output_folder;
            obj.analysis_stage.analysis_outputs = obj.analysis_outputs;
            obj.analysis_stage.extra_analysis_inputs = obj.extra_analysis_inputs;
            obj.analysis_stage.extra_analysis_outputs = obj.extra_analysis_outputs;
            obj.analysis_stage.intermed_result_struct = obj.intermed_result_struct;

            obj.postpro_stage.p = obj.p;
            obj.postpro_stage.b = obj.b;
            obj.postpro_stage.output_folder = obj.output_folder;
            obj.postpro_stage.analysis_outputs = obj.analysis_outputs;
            obj.postpro_stage.postpro_outputs = obj.postpro_outputs;
            obj.postpro_stage.extra_postpro_inputs = obj.extra_postpro_inputs;
            obj.postpro_stage.extra_postpro_outputs = obj.extra_postpro_outputs;
            obj.postpro_stage.fig_array = obj.fig_array;
            obj.postpro_stage.tab_array_display = obj.tab_array_display;
            obj.postpro_stage.tab_array_latex = obj.tab_array_latex;
            obj.postpro_stage.tab_firstrows = obj.tab_firstrows;
            obj.postpro_stage.tab_colspecs = obj.tab_colspecs;
            obj.postpro_stage.end_result_struct = obj.end_result_struct;
        end

        function obj = copy_from_stages(obj)
            obj.p = obj.analysis_stage.p;
            obj.b = obj.analysis_stage.b;
            obj.output_folder = obj.analysis_stage.output_folder;
            obj.intermed_result_struct = obj.analysis_stage.intermed_result_struct;
            obj.analysis_outputs = obj.analysis_stage.analysis_outputs;
            obj.extra_analysis_inputs = obj.analysis_stage.extra_analysis_inputs;
            obj.extra_analysis_outputs = obj.analysis_stage.extra_analysis_outputs;

            obj.fig_array = obj.postpro_stage.fig_array;
            obj.tab_array_display = obj.postpro_stage.tab_array_display;
            obj.tab_array_latex = obj.postpro_stage.tab_array_latex;
            obj.tab_firstrows = obj.postpro_stage.tab_firstrows;
            obj.tab_colspecs = obj.postpro_stage.tab_colspecs;
            obj.end_result_struct = obj.postpro_stage.end_result_struct;
            obj.postpro_outputs = obj.postpro_stage.postpro_outputs;
            obj.extra_postpro_inputs = obj.postpro_stage.extra_postpro_inputs;
            obj.extra_postpro_outputs = obj.postpro_stage.extra_postpro_outputs;
        end

        function obj = run_analysis(obj)
            obj.apply_to_stages();
            obj.analysis_stage = obj.analysis_stage.run_analysis();
            obj.copy_from_stages();
        end

        function obj = run_post_process(obj)
            obj.apply_to_stages();
            obj.postpro_stage = obj.postpro_stage.run_post_process(obj.analysis_stage.intermed_result_struct);
            obj.copy_from_stages();
        end

        function obj = load_intermed_results(obj)
            obj.apply_to_stages();
            obj.analysis_stage = obj.analysis_stage.load_intermed_results();
            obj.copy_from_stages();
        end

        function obj = run_all_from_analysis(obj)
            obj = obj.run_analysis();
            obj = obj.run_post_process();
        end

        function obj = run_all_from_load(obj)
            obj = obj.load_intermed_results();
            obj = obj.run_post_process();
        end

        function obj = run_analysis_from_load_if_possible(obj)
            if isfile([obj.output_folder filesep 'intermed.mat'])
                obj = obj.load_intermed_results();
            else
                obj = obj.run_analysis();
            end
        end

        function obj = run_all_from_load_if_possible(obj)
            if obj.postpro_outputs_exist()
                obj = obj.load_postpro_results();
            else
                obj = obj.run_analysis_from_load_if_possible();
                obj = obj.run_post_process();
            end
        end

        function tf = postpro_outputs_exist(obj)
            obj.apply_to_stages();
            tf = obj.postpro_stage.postpro_outputs_exist();
        end

        function obj = load_postpro_results(obj)
            obj = obj.load_intermed_results();
            obj.apply_to_stages();
            obj.postpro_stage = obj.postpro_stage.load_postpro_results();
            obj.copy_from_stages();
        end

        function stages = write_calkit_stage(obj)
            obj.apply_to_stages();
            analysis_stage = ['analysis-' class(obj) ':' newline ...
                              '  kind: matlab-command' newline ...
                              '  environment: _system' newline ...
                              '  command: add_mdocean_path(); obj=' class(obj) '; obj.run_analysis();' newline ...
                              '  inputs: ' newline ...
                              AbstractStageClass.format_inputs_list(obj.extra_analysis_inputs, obj.analysis_dependencies) newline ...
                              '  outputs: ' newline ...
                              AbstractStageClass.format_outputs([obj.analysis_outputs, obj.extra_analysis_outputs]) ];
            postpro_stage = ['postpro-' class(obj) ':' newline ...
                              '  kind: matlab-command' newline ...
                              '  environment: _system' newline ...
                              '  command: add_mdocean_path(); obj=' class(obj) '; obj.run_all_from_load();' newline ...
                              '  inputs: ' newline ...
                              '    - from_stage_outputs: analysis-' class(obj) newline ...
                              AbstractStageClass.format_inputs_list(obj.extra_postpro_inputs, obj.postpro_dependencies) newline ...
                              '  outputs: ' newline ...
                              AbstractStageClass.format_outputs([obj.postpro_outputs, obj.extra_postpro_outputs]) ];
            viz_stage = ['viz-' class(obj) ':' newline ...
                          '  kind: jupyter-notebook' newline ...
                          '  environment: pubs' newline ...
                          '  notebook_path: results/fig_notebooks/' class(obj) '.ipynb' newline ...
                          '  inputs:' newline ...
                          '    - from_stage_outputs: postpro-' class(obj) newline ...
                          '  html_storage: null' newline ...
                          '  executed_ipynb_storage: git'];
            stages = [analysis_stage newline postpro_stage newline viz_stage];
        end
    end
end
