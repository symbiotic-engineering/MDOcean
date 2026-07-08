classdef (Abstract) AbstractStageClass < handle
    %ABSTRACTSTAGECLASS Common utilities for analysis and postprocessing stages.

    properties
        p
        b
        output_folder
        class_name
        class_dependencies_cached
    end

    properties (Dependent)
        class_dependencies
    end

    methods
        function obj = AbstractStageClass(class_name, p, b)
            output_folder_rel = ['results' filesep class_name];
            MDOcean_folder = fileparts(fileparts(fileparts(mfilename('fullpath'))));
            obj.output_folder = fullfile(MDOcean_folder, output_folder_rel);
            if ~isfolder(obj.output_folder)
                mkdir(obj.output_folder)
            end

            obj.class_name = class_name;
            if nargin > 1 && ~isempty(p)
                obj.p = p;
            else
                obj.p = parameters();
            end
            if nargin > 2 && ~isempty(b)
                obj.b = b;
            else
                if isfile('var_bounds.mat')
                    obj.b = load('var_bounds.mat').b;
                else
                    obj.b = var_bounds();
                end
            end
        end

        function val = get.class_dependencies(obj)
            if isempty(obj.class_dependencies_cached)
                obj.class_dependencies_cached = obj.get_dependencies(class(obj));
            end
            val = obj.class_dependencies_cached;
        end
    end

    methods (Static)
        function deps_rel = get_dependencies(fcn_name)
            fcn_deps = matlab.codetools.requiredFilesAndProducts(fcn_name);
            fcn_path = which(fcn_name);
            if isempty(fcn_path)
                deps_abs = fcn_deps;
            else
                deps_abs = [fcn_deps, {fcn_path}];
            end

            path = mfilename('fullpath');
            MDOcean_folder = fileparts(fileparts(fileparts(path)));
            deps_rel = AbstractStageClass.make_rel_path(deps_abs, MDOcean_folder);
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

        function list_str = format_inputs_list(extras, deps)
            combined = [string(extras), string(deps)];
            combined_unique = unique(combined, 'stable');
            lines = "    - " + combined_unique;
            list_str = char(join(lines, newline));
        end

        function list_str = format_outputs(outputs)
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
    end
end
