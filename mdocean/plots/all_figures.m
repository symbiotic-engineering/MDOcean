% Generate all figures used in the paper
function [fig_success,tab_success,...
          fig_output, tab_output,...
          fig_runtime,tab_runtime,...
          num_figs,num_tabs,...
          fig_names,tab_names] = all_figures( which_figs, which_tabs, filename_uuid )

if nargin<3
    filename_uuid = ''; % required argument for anything running gradient_optim in parallel,
                        % since prob2struct needs unique filenames for code generation
end

date = datestr(now,'yyyy-mm-dd_HH.MM.SS');
save_folder = ['../test-results/' date '/'];
mkdir(save_folder)
table_save_fcn = @(tab,filename,varargin) table2latex(tab, [save_folder filename], varargin{:});

num_figs = 59;
num_tabs = 8;

if nargin==0
    % if run without arguments, show all figures and tables
    which_figs = 1:num_figs;
    which_tabs = 1:num_tabs;
end

fig_success = cell([1,length(which_figs)]);
tab_success = cell([1,length(which_tabs)]);

fig_output = gobjects(1, length(which_figs));
tab_output(1, 1:length(which_tabs)) = {table()};

fig_runtime = NaN([1,length(which_figs)]);
tab_runtime = NaN([1,length(which_tabs)]);
%%



%% Auto-generate figure and table names from list above
fig_nums = cellstr(num2str((1:num_figs).'));
tab_idxs = ~cellfun(@isempty,tabs_in_paper);
tab_nums = cellstr(num2str(find(tab_idxs).'));

fig_names = strcat("Fig. ", fig_nums, ": ", figs_in_paper.');
tab_names = strcat("Tab. ", tab_nums, ": ", tabs_in_paper(tab_idxs).');

%% Initialize structures to store generated figures and tables
generated_figs = struct();
generated_tabs = struct();

% Create analysis class instances cache
analysis_instances = containers.Map();

% Loop for figures
for i = 1:length(which_figs)
    try
        % Get the class name for the figure
        fig_number = which_figs(i);
        fig_name = split(figs_in_paper{fig_number}, '.');
        class_name = fig_name{1};
        fig_desc = fig_name{2};
        
        % Check if the figure has already been generated
        t = tic;
        if ~isfield(generated_figs, class_name)
            % Create or retrieve analysis instance
            if ~isKey(analysis_instances, class_name)
                % Create new instance and set filename_uuid and table_save_fcn
                analysis_obj = eval([class_name '()']);
                analysis_obj.b.filename_uuid = filename_uuid;
                analysis_obj.b.table_save_fcn = table_save_fcn;
                analysis_instances(class_name) = analysis_obj;
            else
                analysis_obj = analysis_instances(class_name);
            end
            
            % Run the analysis if not already done
            if isempty(analysis_obj.fig_array)
                analysis_obj = analysis_obj.run_analysis();
                analysis_obj = analysis_obj.run_post_process();
                analysis_instances(class_name) = analysis_obj;
            end
            
            % Extract figures and tables from the analysis object
            figs = struct();
            tabs = struct();
            for j = 1:length(analysis_obj.fig_names)
                if j <= length(analysis_obj.fig_array)
                    figs.(analysis_obj.fig_names{j}) = analysis_obj.fig_array(j);
                end
            end
            for j = 1:length(analysis_obj.tab_names)
                if j <= length(analysis_obj.tab_array_display)
                    tabs.(analysis_obj.tab_names{j}) = analysis_obj.tab_array_display{j};
                end
            end
            
            % Store the generated figure in the generated_figs structure
            generated_figs.(class_name) = figs;
            
            % If the function generates a table, store it in the generated_tabs structure
            if ~isempty(tabs)
                generated_tabs.(class_name) = tabs;
            end
        end
        
        % Store the generated figure in fig_output
        fig_output(i) = generated_figs.(class_name).(fig_desc);
        fig_runtime(i) = toc(t);

    catch err
        fig_success{i} = err;  % Store error for the figure
    end
end

% Loop for tables
for i = 1:length(which_tabs)
    try
        tab_number = which_tabs(i);
        % Get the class name for the table
        tab_name = split(tabs_in_paper{tab_number}, '.');
        class_name = tab_name{1};
        tab_desc = tab_name{2};
        
        % Check if the table has already been generated
        t = tic;
        if ~isfield(generated_tabs, class_name)
            % Create or retrieve analysis instance
            if ~isKey(analysis_instances, class_name)
                % Create new instance and set filename_uuid and table_save_fcn
                eval(['analysis_obj = ' class_name '();']);
                analysis_obj.b.filename_uuid = filename_uuid;
                analysis_obj.b.table_save_fcn = table_save_fcn;
                analysis_instances(class_name) = analysis_obj;
            else
                analysis_obj = analysis_instances(class_name);
            end
            
            % Run the analysis if not already done
            if isempty(analysis_obj.tab_array_display)
                analysis_obj = analysis_obj.run_analysis();
                analysis_obj = analysis_obj.run_post_process();
                analysis_instances(class_name) = analysis_obj;
            end
            
            % Extract figures and tables from the analysis object
            figs = struct();
            tabs = struct();
            for j = 1:length(analysis_obj.fig_names)
                if j <= length(analysis_obj.fig_array)
                    figs.(analysis_obj.fig_names{j}) = analysis_obj.fig_array(j);
                end
            end
            for j = 1:length(analysis_obj.tab_names)
                if j <= length(analysis_obj.tab_array_display)
                    tabs.(analysis_obj.tab_names{j}) = analysis_obj.tab_array_display{j};
                end
            end
            
            % Store the generated table in the generated_tabs structure
            generated_tabs.(class_name) = tabs;
            
            % If the function generates a figure, store it in the generated_figs structure
            if ~isempty(figs)
                generated_figs.(class_name) = figs;
            end
        end
        
        % Store the generated table in tab_output
        tab_output{i} = generated_tabs.(class_name).(tab_desc);  % Store the table
        tab_runtime(i) = toc(t);
        display(tab_output{i})
        
    catch err
        tab_success{i} = err;  % Store error for the table
    end
end

end
