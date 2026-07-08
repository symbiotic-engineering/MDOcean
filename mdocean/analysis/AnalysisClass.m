classdef AnalysisClass < AbstractStageClass
    %ANALYSISCLASS Analysis-stage execution and intermediate caching.

    properties
        analysis_fcn_handle
        intermed_result_struct
        analysis_outputs
        extra_analysis_inputs = {'./mdocean/inputs/validation/RM3-CBS.xlsx'}
        extra_analysis_outputs = {}
    end

    properties (Dependent)
        analysis_dependencies
    end

    methods
        function obj = AnalysisClass(class_name, analysis_fcn_handle, p, b)
            obj@AbstractStageClass(class_name, p, b)
            obj.analysis_fcn_handle = analysis_fcn_handle;
            output_folder_rel = ['results' filesep class_name];
            obj.analysis_outputs = {[output_folder_rel, filesep, 'intermed.mat']};
        end

        function val = get.analysis_dependencies(obj)
            analysis_deps = obj.get_dependencies(['@' obj.class_name filesep 'analysis_fcn']);
            stage_deps = {'./mdocean/analysis/AbstractStageClass.m', './mdocean/analysis/AnalysisClass.m'};
            all_deps = [stage_deps, analysis_deps];
            sorted = sort(unique(all_deps));
            to_remove = 'OpenFLASH';
            val = sorted(~contains(sorted, to_remove));
        end

        function obj = run_analysis(obj)
            original_folder = pwd;
            mdocean_folder = fileparts(fileparts(mfilename('fullpath')));
            cleanup_obj = onCleanup(@() cd(original_folder)); %#ok<NASGU>
            cd(mdocean_folder);
            t = tic;
            intermed_result_struct = obj.analysis_fcn_handle(obj.p, obj.b);
            intermed_result_struct.analysis_time = toc(t);
            obj.intermed_result_struct = intermed_result_struct;
            obj.save_intermed_results();
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
    end
end
