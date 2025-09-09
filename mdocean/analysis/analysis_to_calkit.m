analysis_list = dir('analysis/*.m');
analysis_list = {analysis_list.name};
analysis_list = analysis_list(~contains(analysis_list, ["GenericAnalysis.m","analysis_to_calkit.m"]));

stages = cell(length(analysis_list),1);

p = parameters();
b = var_bounds();

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