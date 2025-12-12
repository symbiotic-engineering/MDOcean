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
                          strcat(obj.tab_names,'.tex') ...
                          'end.mat']);
        end
        function val = get.class_dependencies(obj)
            val = obj.get_dependencies(class(obj));
        end
        function val = get.analysis_dependencies(obj)
            all_deps = obj.class_dependencies;
            to_remove = {['@' class(obj) filesep 'post_process_fcn'], 'OpenFLASH'};
            val = all_deps(~contains(all_deps, to_remove));
        end
        function val = get.postpro_dependencies(obj)
            all_deps = obj.class_dependencies;
            to_remove = {['@' class(obj) filesep 'analysis_fcn'], 'OpenFLASH'};
            val = all_deps(~contains(all_deps, to_remove));
        end
        function obj = run_analysis(obj)
            cd('mdocean');
            obj.intermed_result_struct = obj.analysis_fcn(obj.p, obj.b);
            cd('..');
            obj.save_intermed_results();
        end

        function obj = run_post_process(obj)
            cd('mdocean');
            [obj.fig_array,...
                obj.tab_array_display,...
                obj.tab_array_latex,...
                obj.end_result_struct,...
                obj.tab_firstrows,...
                obj.tab_colspecs] = obj.post_process_fcn(obj.intermed_result_struct);
            cd('..');

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
            cell2filelist = @(c) char(join(strcat("    - ",c),newline));

            analysis_stage = ['analysis-' class(obj) ':' newline ...
                              '  kind: matlab-command' newline ...
                              '  environment: _system' newline ...
                              '  command: add_mdocean_path(); obj=' class(obj) '; obj.run_analysis();' newline ...
                              '  inputs: ' newline ...
                              cell2filelist(obj.analysis_dependencies.') newline ...
                              '  outputs: ' newline ...
                              cell2filelist(obj.analysis_outputs.') ];
            postpro_stage = ['postpro-' class(obj) ':' newline ...
                              '  kind: matlab-command' newline ...
                              '  environment: _system' newline ...
                              '  command: add_mdocean_path(); obj=' class(obj) '; obj.run_all_from_load();' newline ...
                              '  inputs: ' newline ...
                              cell2filelist(obj.postpro_dependencies.') newline ...
                              '  outputs: ' newline ...
                              cell2filelist(obj.postpro_outputs.') ];
            stages = [analysis_stage newline postpro_stage];
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

            valid_mask = ishghandle(figs_flat) & isprop(figs_flat,'Type') & strcmp({figs_flat.Type},'figure');

            figs_flat(~valid_mask) = gobjects(1);

            % assign to output, taking into account if figs_flat is smaller or larger than desired
            figs_out = gobjects(1,length(obj.fig_names));
            len = min(length(figs_flat), length(figs_out));
            figs_out(1:len) = figs_flat(1:len);
        end
    end
    methods (Static)
        function deps_rel = get_dependencies(fcn_name)
            deps_abs = matlab.codetools.requiredFilesAndProducts(fcn_name);
            path = mfilename('fullpath');
            s = split(which(path), filesep);
            MDOcean_folder = strjoin(s(1:end-3), filesep);
            base = MDOcean_folder
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