classdef PostproTableClass < PostproClass
    %POSTPROTABLECLASS Postprocessing class with LaTeX table outputs.

    methods
        function obj = PostproTableClass(class_name, post_process_fcn_handle, fig_names, tab_names, p, b)
            obj@PostproClass(class_name, post_process_fcn_handle, fig_names, tab_names, p, b)
        end

        function outputs = build_postpro_outputs(obj, output_folder_rel)
            outputs = strcat(output_folder_rel,...
                filesep, [strcat(obj.fig_names,'.pdf') ...
                          strcat(obj.fig_names,'.fig') ...
                          strcat(obj.tab_names,'.tex') ...
                          'end.mat' ...
                          'end.json']);
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
    end
end
