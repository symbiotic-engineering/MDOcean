classdef (Abstract) GenericAnalysis
    %GENERICANALYSIS Abstract class for MDOcean analyses
    %   Handles saving and loading of intermediate results, figures, and tables.

    properties
        p
        b
        fig_array
        tab_array_display
        tab_array_latex
        tab_firstrows
        tab_colspecs
        intermed_result_struct % all results needed to debug and make plots
        end_result_struct % only summary end results for latex output or tests
        output_folder
        analysis_outputs
        postpro_outputs
        extra_analysis_inputs = {'./mdocean/inputs/validation/RM3-CBS.xlsx'} % non-code data inputs for calkit analysis stage
        extra_analysis_outputs = {} % non-code extra outputs beyond intermed.mat for calkit analysis stage
        extra_postpro_inputs = {'./mdocean/inputs/validation/RM3-CBS.xlsx'} % non-code data inputs for calkit postpro stage
        extra_postpro_outputs = {} % non-code extra outputs for calkit postpro stage
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
        function obj = GenericAnalysis(p,b)
            obj.output_folder = ['results' filesep class(obj)];
            if ~isfolder(obj.output_folder)
                mkdir(obj.output_folder)
            end
            if nargin > 0
                obj.p = p;
            else
                obj.p = parameters();
            end
            if nargin > 1
                obj.b = b;
            else
                obj.b = var_bounds();
            end
            
            obj.analysis_outputs = {[obj.output_folder, filesep, 'intermed.mat']};

            obj.postpro_outputs  = strcat(obj.output_folder,...
                filesep, [strcat(obj.fig_names,'.pdf') ...
                          strcat(obj.fig_names,'.fig') ...
                          strcat(obj.tab_names,'.tex') ...
                          'end.mat' ...
                          'end.json']);
        end
        function val = get.class_dependencies(obj)
            parent = superclasses(obj);
            val = obj.get_dependencies(parent{:});
        end
        function val = get.analysis_dependencies(obj)
            analysis_deps = obj.get_dependencies(['@' class(obj) filesep 'analysis_fcn']);
            all_deps = [obj.class_dependencies, analysis_deps];
            sorted = sort(unique(all_deps));
            to_remove = 'OpenFLASH';
            val = sorted(~contains(sorted, to_remove));
        end
        function val = get.postpro_dependencies(obj)
            postpro_deps = obj.get_dependencies(['@' class(obj) filesep 'post_process_fcn']);
            all_deps = [obj.class_dependencies, postpro_deps, obj.analysis_outputs];
            sorted = sort(unique(all_deps));
            to_remove = 'OpenFLASH';
            val = sorted(~contains(sorted, to_remove));
        end
        function obj = run_analysis(obj)
            cd('mdocean');
            t = tic;
            intermed_result_struct = obj.analysis_fcn(obj.p, obj.b);
            intermed_result_struct.analysis_time = toc(t);
            obj.intermed_result_struct = intermed_result_struct;
            cd('..');
            obj.save_intermed_results();
        end

        function obj = run_post_process(obj)
            cd('mdocean');
            t = tic;
            [obj.fig_array,...
                obj.tab_array_display,...
                obj.tab_array_latex,...
                end_result_struct,...
                obj.tab_firstrows,...
                obj.tab_colspecs] = obj.post_process_fcn(obj.intermed_result_struct);
            cd('..');
            end_result_struct.postpro_time = toc(t);
            end_result_struct.analysis_time = obj.intermed_result_struct.analysis_time;
            obj.end_result_struct = end_result_struct;

            obj.fig_array = obj.validate_figs(obj.fig_array);

            obj.save_figs();
            obj.save_tabs();
            obj.save_end_results();
        end

        function obj = load_intermed_results(obj)
            fname = [obj.output_folder filesep 'intermed.mat'];
            if ~isfile(fname)
                error('No intermediate results file found. Run analysis first.')
            end
            s = load(fname);
            d = dir([obj.output_folder filesep 'intermed_*.fig']);
            fnames = {d.name};
            tokens = regexp(fnames, 'intermed_(.+?)_(\d+).fig', 'tokens');
            [var_names,~,idx] = unique(cellfun(@(t) t{1}{1}, tokens, 'UniformOutput', false));
            for i=1:length(var_names)
                var_name = var_names{i};
                fig_idxs = extractBetween(fnames(idx==i), ['intermed_' var_name '_'], '.fig');
                fig_idxs = str2double(fig_idxs);
                fig_vec = gobjects(1,length(fig_idxs));
                for j=1:length(fig_idxs)
                    fig_name = ['intermed_' var_name '_' num2str(fig_idxs(j)) '.fig'];
                    fig_path = [obj.output_folder filesep fig_name];
                    try
                        fig_handle = openfig(fig_path, 'invisible');
                        if isfield(fig_handle.UserData, 'Position')
                            fig_handle.Position(3:4) = fig_handle.UserData.Position;
                        end
                        fig_vec(j) = fig_handle;
                    catch
                        warning(['Could not load figure: ' fig_path])
                    end
                end
                s.(var_name) = fig_vec;
            end
            obj.intermed_result_struct = s;
        end

        function save_intermed_results(obj)
            s = obj.intermed_result_struct;
            fn = fieldnames(s);
            for i=1:length(fn)
                if any(isgraphics(s.(fn{i})) & isa(s.(fn{i}),'matlab.ui.Figure'))
                    figs = s.(fn{i});
                        for j=1:length(figs)
                            fig = figs(j);
                            fig = check_fig_size(fig);
                            savefig(fig, [obj.output_folder filesep 'intermed_', fn{i}, '_', num2str(j), '.fig'])
                        end
                    s = rmfield(s, fn{i});
                end
            end
            save([obj.output_folder filesep 'intermed'], '-struct', 's')
        end

        function save_end_results(obj)
            s =  obj.end_result_struct;
            fname = [obj.output_folder filesep 'end'];
            save(fname,'-struct', 's')
            json = jsonencode(s,PrettyPrint=true);
            writelines(json, [fname '.json']);
        end

        function save_figs(obj)
            for i=1:length(obj.fig_names)
                fig = obj.fig_array(i);
                save_fig_with_diagnostic(fig, obj.fig_names{i}, obj.output_folder)
            end
        end
        
        function save_tabs(obj)
            for i=1:length(obj.tab_names)
                tab = obj.tab_array_latex{i};
                fname = [obj.output_folder filesep obj.tab_names{i}];
                colspec = obj.tab_colspecs{i};
                firstrow = obj.tab_firstrows{i};
                table2latex(tab,fname,colspec,firstrow)
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

        function obj = run_analysis_from_load_if_possible(obj)
            if isempty(obj.intermed_result_struct)
                obj = obj.run_analysis();
            else
                obj = obj.load_intermed_results();
            end
        end

        function obj = run_all_from_load_if_possible(obj)
            obj = obj.run_analysis_from_load_if_possible();
            obj = obj.run_post_process();
        end

        function stages = write_calkit_stage(obj)
            analysis_stage = ['analysis-' class(obj) ':' newline ...
                              '  kind: matlab-command' newline ...
                              '  environment: _system' newline ...
                              '  command: add_mdocean_path(); obj=' class(obj) '; obj.run_analysis();' newline ...
                              '  inputs: ' newline ...
                              obj.format_inputs_list(obj.extra_analysis_inputs, obj.analysis_dependencies) newline ...
                              '  outputs: ' newline ...
                              obj.format_outputs([obj.analysis_outputs, obj.extra_analysis_outputs]) ];
            postpro_stage = ['postpro-' class(obj) ':' newline ...
                              '  kind: matlab-command' newline ...
                              '  environment: _system' newline ...
                              '  command: add_mdocean_path(); obj=' class(obj) '; obj.run_all_from_load();' newline ...
                              '  inputs: ' newline ...
                              '    - from_stage_outputs: analysis-' class(obj) newline ...
                              obj.format_inputs_list(obj.extra_postpro_inputs, obj.postpro_dependencies) newline ...
                              '  outputs: ' newline ...
                              obj.format_outputs([obj.postpro_outputs, obj.extra_postpro_outputs]) ];
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

        function list_str = format_inputs_list(~, extras, deps)
            %FORMAT_INPUTS_LIST Format inputs as YAML list: extras first, then deps, de-duplicated.
            combined = [string(extras), string(deps)];
            combined_unique = unique(combined, 'stable');
            lines = "    - " + combined_unique;
            list_str = char(join(lines, newline));
        end

        function list_str = format_outputs(~, outputs)
            %FORMAT_OUTPUTS Format a cell array of output paths as YAML list entries.
            %   .tex and .json outputs include 'storage: git'; others are plain paths.
            outputs = string(outputs);
            if isempty(outputs)
                list_str = '';
                return;
            end
            [~, ~, exts] = fileparts(outputs);
            is_git = ismember(exts, [".tex", ".json"]);
            lines = strings(size(outputs));
            lines(~is_git) = "    - " + outputs(~is_git);
            lines(is_git) = "    - path: " + outputs(is_git) + newline + "      storage: git";
            list_str = char(join(lines, newline));
        end

        function figs_out = validate_figs(obj, figs_in)
            %VALIDATE_FIGS Return only valid MATLAB figure handles from input
            % Accepts graphics handles, cell arrays of handles, or empty.
            if isempty(figs_in)
                figs_out = gobjects(0);
                return;
            end

            % Unwrap cell arrays
            if iscell(figs_in)
                try
                    figs_in = [figs_in{:}];
                catch
                    figs_out = gobjects(0);
                    return;
                end
            end

            % Flatten
            figs_flat = figs_in(:).';

            % check if figure
            has_type = isprop(figs_flat,'Type');
            is_fig = false(size(figs_flat));
            is_fig(has_type) = strcmp({figs_flat(has_type).Type},'figure');
            valid_mask = ishghandle(figs_flat) &  is_fig;

            figs_flat(~valid_mask) = gobjects(1);

            % assign to output, taking into account if figs_flat is smaller or larger than desired
            figs_out = gobjects(1,length(obj.fig_names));
            len = min(length(figs_flat), length(figs_out));
            figs_out(1:len) = figs_flat(1:len);
        end
    end
    methods (Static)
        function deps_rel = get_dependencies(fcn_name)
            fcn_deps = matlab.codetools.requiredFilesAndProducts(fcn_name);
            deps_abs = [fcn_deps, which(fcn_name)];
            path = mfilename('fullpath');
            s = split(which(path), filesep);
            MDOcean_folder = strjoin(s(1:end-3), filesep);
            base = MDOcean_folder;
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
