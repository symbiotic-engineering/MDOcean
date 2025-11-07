
p = parameters();
b = var_bounds();

analysis_list = get_all_analyses();
stages = cell(length(analysis_list),1);

for i = 1:length(analysis_list)
    [~, class_name] = fileparts(analysis_list{i});

    try
        analysis_obj = feval(class_name,p,b);
        stage = analysis_obj.write_calkit_stage();
    catch ME
        warning("Failed to write calkit stage for %s: %s", class_name, ME.message);
        stage = '';
    end
    stages{i} = stage;
end

disp(stages)
stages_combined = strjoin(stages,newline);
disp(stages_combined)

writelines(stages_combined,'calkit_stages.txt') % then manually copy-paste into calkit.yaml
