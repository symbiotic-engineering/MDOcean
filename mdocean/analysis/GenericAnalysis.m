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
        function obj = GenericAnalysis(p,b)
            obj.output_folder = ['results' filesep class(obj)];
            if ~isfolder(obj.output_folder)
                mkdir(obj.output_folder)
            end
            if nargin > 0
                obj.p = p;
            end
            if nargin > 1
                obj.b = b;
            end
            obj.analysis_dependencies = obj.get_dependencies([class(obj) '.analysis_fcn()']);
            obj.analysis_outputs = {[obj.output_folder, filesep, 'intermed.mat']};
            obj.postpro_dependencies = [obj.get_dependencies([class(obj) '.post_process_fcn()']),...
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
            fname = [obj.output_folder filesep 'intermed.mat'];
            if ~isfile(fname)
                error('No intermediate results file found. Run analysis first.')
            end
            obj.intermed_result_struct = load(fname);
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

        function obj = run_all_from_analysis(obj)
            obj = obj.run_analysis();
            obj = obj.run_post_process();
        end

        function obj = run_all_from_load(obj)
            obj = obj.load_intermed_results();
            obj = obj.run_post_process();
        end

        function stages = write_calkit_stage(obj)
            cell2filelist = @(c) char(join(strcat("    - ",c),newline));

            analysis_stage = ['analysis-' class(obj) ':' newline ...
                              '  kind: matlab-command' newline ...
                              '  environment: matlab' newline ...
                              '  command: add_mdocean_path(); obj=' class(obj) '; obj.run_analysis();' newline ...
                              '  inputs: ' newline ...
                              cell2filelist(obj.analysis_dependencies.') newline ...
                              '  outputs: ' newline ...
                              cell2filelist(obj.analysis_outputs.') ];
            postpro_stage = ['postpro-' class(obj) ':' newline ...
                              '  kind: matlab-command' newline ...
                              '  environment: matlab' newline ...
                              '  command: add_mdocean_path(); obj=' class(obj) '; obj.run_all_from_load();' newline ...
                              '  inputs: ' newline ...
                              cell2filelist(obj.postpro_dependencies.') newline ...
                              '  outputs: ' newline ...
                              cell2filelist(obj.postpro_outputs.') ];
            stages = [analysis_stage newline postpro_stage];
        end
    end
    methods (Static)
        function deps_rel = get_dependencies(fcn_name)
            deps_abs = matlab.codetools.requiredFilesAndProducts(fcn_name);
            base = what('..').path;
            deps_rel = GenericAnalysis.make_rel_path(deps_abs, base);
        end
        
        function rel_paths = make_rel_path(abs_paths, base)
            rel_paths = cell(size(abs_paths));
            for i = 1:numel(abs_paths)
                [path_dir, file_name, file_ext] = fileparts(abs_paths{i});

                if startsWith(path_dir, base)
                    rel_path = strrep(path_dir, base, '.');
                    rel_path = fullfile(rel_path, strcat(file_name, file_ext));
                else
                    error('Dependency outside of base folder')
                end
                rel_paths{i} = rel_path;
            end
        end
    end
end