classdef PostproClass < AbstractStageClass
    %POSTPROCLASS Postprocessing-stage execution and final artifact caching.

    properties
        post_process_fcn_handle
        fig_names
        tab_names
        fig_array
        tab_array_display
        tab_array_latex
        tab_firstrows
        tab_colspecs
        end_result_struct
        analysis_outputs
        postpro_outputs
        extra_postpro_inputs = {'./mdocean/inputs/validation/RM3-CBS.xlsx'}
        extra_postpro_outputs = {}
    end

    properties (Dependent)
        postpro_dependencies
    end

    methods
        function obj = PostproClass(class_name, post_process_fcn_handle, fig_names, tab_names, p, b)
            obj@AbstractStageClass(class_name, p, b)
            obj.post_process_fcn_handle = post_process_fcn_handle;
            obj.fig_names = fig_names;
            obj.tab_names = tab_names;
            output_folder_rel = ['results' filesep class_name];
            obj.analysis_outputs = {[output_folder_rel, filesep, 'intermed.mat']};
            obj.postpro_outputs = obj.build_postpro_outputs(output_folder_rel);
        end

        function val = get.postpro_dependencies(obj)
            postpro_deps = obj.get_dependencies(['@' obj.class_name filesep 'post_process_fcn']);
            if isa(obj, 'PostproTableClass')
                stage_class_dep = './mdocean/analysis/PostproTableClass.m';
            else
                stage_class_dep = './mdocean/analysis/PostproClass.m';
            end
            stage_deps = {'./mdocean/analysis/AbstractStageClass.m', stage_class_dep};
            all_deps = [stage_deps, postpro_deps, obj.analysis_outputs];
            sorted = sort(unique(all_deps));
            to_remove = 'OpenFLASH';
            val = sorted(~contains(sorted, to_remove));
        end

        function outputs = build_postpro_outputs(obj, output_folder_rel)
            outputs = strcat(output_folder_rel,...
                filesep, [strcat(obj.fig_names,'.pdf') ...
                          strcat(obj.fig_names,'.fig') ...
                          'end.mat' ...
                          'end.json']);
        end

        function obj = run_post_process(obj, intermed_result_struct)
            original_folder = pwd;
            mdocean_folder = fileparts(fileparts(mfilename('fullpath')));
            cleanup_obj = onCleanup(@() cd(original_folder)); %#ok<NASGU>
            cd(mdocean_folder);
            t = tic;
            [obj.fig_array,...
                obj.tab_array_display,...
                obj.tab_array_latex,...
                end_result_struct,...
                obj.tab_firstrows,...
                obj.tab_colspecs] = obj.post_process_fcn_handle(intermed_result_struct);
            end_result_struct.postpro_time = toc(t);
            if isfield(intermed_result_struct, 'analysis_time')
                end_result_struct.analysis_time = intermed_result_struct.analysis_time;
            end
            obj.end_result_struct = end_result_struct;

            obj.fig_array = obj.validate_figs(obj.fig_array);

            obj.save_figs();
            obj.save_tabs();
            obj.save_end_results();
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

        function save_tabs(~)
            % No table output for generic postprocessing class.
        end

        function tf_2 = postpro_outputs_exist(obj)
            fig_paths = {};
            tab_paths = {};
            if ~isempty(obj.fig_names)
                fig_paths = fullfile(obj.output_folder, strcat(obj.fig_names, '.pdf'));
            end
            if ~isempty(obj.tab_names)
                tab_paths = fullfile(obj.output_folder, strcat(obj.tab_names, '.tex'));
            end
            required = [{fullfile(obj.output_folder, 'intermed.mat')}, ...
                        {fullfile(obj.output_folder, 'end.mat')}, ...
                        fig_paths, ...
                        tab_paths];
            tf = all(cellfun(@isfile, required));
            tf_2 = tf && isempty(obj.tab_names);
        end

        function obj = load_postpro_results(obj)
            s = load(fullfile(obj.output_folder, 'end.mat'));
            obj.end_result_struct = s;

            num_figs = length(obj.fig_names);
            fig_array = gobjects(1, num_figs);
            for i = 1:num_figs
                fig_path = fullfile(obj.output_folder, [obj.fig_names{i} '.fig']);
                try
                    fig_handle = openfig(fig_path, 'invisible');
                    if ishghandle(fig_handle) ...
                            && isstruct(fig_handle.UserData) ...
                            && isfield(fig_handle.UserData, 'Position') ...
                            && numel(fig_handle.UserData.Position) == 2
                        fig_handle.Position(3:4) = fig_handle.UserData.Position;
                    end
                    fig_array(i) = fig_handle;
                catch err
                    warning('Failed to load cached figure from %s: %s.', fig_path, err.message)
                end
            end
            obj.fig_array = fig_array;
        end

        function figs_out = validate_figs(obj, figs_in)
            if isempty(figs_in)
                figs_out = gobjects(0);
                return;
            end

            if iscell(figs_in)
                try
                    figs_in = [figs_in{:}];
                catch
                    figs_out = gobjects(0);
                    return;
                end
            end

            figs_flat = figs_in(:).';

            has_type = isprop(figs_flat,'Type');
            is_fig = false(size(figs_flat));
            is_fig(has_type) = strcmp({figs_flat(has_type).Type},'figure');
            valid_mask = ishghandle(figs_flat) &  is_fig;

            figs_flat(~valid_mask) = gobjects(1);

            figs_out = gobjects(1,length(obj.fig_names));
            len = min(length(figs_flat), length(figs_out));
            figs_out(1:len) = figs_flat(1:len);
        end
    end
end
