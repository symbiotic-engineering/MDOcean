function [f1, f2, f3, f4] = module_runtime_compare(p,b)
    
    X = [b.X_noms; 1];
    num_outputs = 3;

    % run sim with timing
    p.use_multibody = true;
    profile clear
    profile on
    simulation(X,p);
    profile_multibody = profile('info');
    t_multibody_fullsim_timeit = timeit(@()simulation(X,p),num_outputs);

    p.use_multibody = false;
    profile clear
    profile on
    simulation(X,p);
    profile_singlebody = profile('info');
    t_singlebody_fullsim_timeit = timeit(@()simulation(X,p),num_outputs);

    % extract full sim timing and multiplier
    t_multibody_fullsim_profile  = extract_runtime(profile_multibody, 'simulation');
    t_singlebody_fullsim_profile = extract_runtime(profile_singlebody,'simulation');
   
    time_multiplier_multibody  = t_multibody_fullsim_timeit / t_multibody_fullsim_profile;
    time_multiplier_singlebody = t_singlebody_fullsim_timeit / t_singlebody_fullsim_profile;

    % just dynamics
    t_multibody  = time_multiplier_multibody  * extract_runtime(profile_multibody, 'get_response_drag');
    t_singlebody = time_multiplier_singlebody * extract_runtime(profile_singlebody,'get_response_drag');

    gcp;
    time = tic;
    pause(1)
%     runRM3Parallel
    t_wecsim = 447; %toc(time); %447

    f1 = figure;
    cats = {'MDOcean Dynamics 1-DOF','MDOcean Dynamics 2-DOF','WecSim (Parallelized)'};
    times = [t_singlebody t_multibody t_wecsim];
    cats = reordercats(categorical(cats),cats);
    h = bar(cats,times);
    if ~isMATLABReleaseOlderThan('R2024b')
        h(1).Labels = h(1).YData;
    end
    title('Dynamics Runtime Comparison')
    subtitle('for all sea states')
    ylabel('Runtime (s)')
    h.Parent.YScale = 'log';
    h.Parent.YGrid  = 'on';
    improvePlot
    h.Parent.YMinorGrid = 'off';
    f1.Position = [100 100 600 640];
    exponents = log10(times);
    ylim(10.^[min(floor(exponents)),max(ceil(exponents))]) % one order of magnitude above max and below min

    % hydro
    % determine line numbers
    meem_filename = 'A_b_c_matrix_N10_M10_K10_heaving_outer';
    bessel_lines     = findStringInFile(meem_filename, 'bessel');
    unpack_lines     = findStringInFile(meem_filename, 'ct{:}');
    dispersion_lines = findStringInFile('run_MEEM','m_k_h_deg = fzero(eqn, bounds)');
    lin_solve_lines  = findStringInFile('run_MEEM','A_num\b_num');

    % ensure no double counting lines from same file
    assert(~any(ismember(bessel_lines,unpack_lines)))
    assert(~any(ismember(dispersion_lines,lin_solve_lines)))

    % extract times for line numbers
    t_hydro            = time_multiplier_singlebody * extract_runtime(profile_singlebody,'get_dynamic_coeffs');
    t_hydro_bessels    = time_multiplier_singlebody * extract_runtime(profile_singlebody,'ft_1',bessel_lines);
    t_hydro_unpack     = time_multiplier_singlebody * extract_runtime(profile_singlebody,'ft_1',unpack_lines);
    t_hydro_dispersion = time_multiplier_singlebody * extract_runtime(profile_singlebody,'compute_eigen_hydro_coeffs',dispersion_lines);
    t_hydro_lin_solve  = time_multiplier_singlebody * extract_runtime(profile_singlebody,'compute_eigen_hydro_coeffs',lin_solve_lines);

    t_hydro_other = t_hydro - (t_hydro_bessels + t_hydro_unpack + t_hydro_dispersion + t_hydro_lin_solve);
    times_meem = [t_hydro_bessels t_hydro_unpack t_hydro_dispersion t_hydro_lin_solve t_hydro_other];

    t_capytaine = 0.323;

    sections_hydro = {'Bessel Functions','Unpacking Variables','Dispersion Relation','Linear Solve','Other'};
    num_freq_hydro = length([p.T p.T_struct]);
    times_hydro = [times_meem/num_freq_hydro 0; zeros(1,length(times_meem)) t_capytaine];
    cats_hydro = {'MDOcean MEEM','Capytaine'};
    cats_hydro = reordercats(categorical(cats_hydro),cats_hydro);

    % hydro figure, not broken down, log scale
    f2 = figure;
    h = barh(cats_hydro,sum(times_hydro,2));
    title('Hydrodynamics Runtime Comparison')
    subtitle('for a single frequency')
    ylabel('Runtime (s)')
    if ~isMATLABReleaseOlderThan('R2024b')
        h(1).Labels = h(1).YData;
    end
    ax = h(1).Parent;
    split_labels_by_space(ax)
    h.Parent.XScale = 'log';
    h.Parent.YGrid  = 'on';
    improvePlot
    h.Parent.YMinorGrid = 'off';
    f2.Position = [100 100 1410 600];

    
    % hydro figure, broken down by category, not log scale
    f3 = figure;
    h = barh(cats_hydro,times_hydro,'stacked');
    legend(sections_hydro)
    title('Hydrodynamics Runtime Comparison')
    subtitle('for a single frequency')
    ylabel('Runtime (s)')
    if ~isMATLABReleaseOlderThan('R2024b')
        h(1).Labels = h(1).YData;
    end
    ax = h(1).Parent;
    split_labels_by_space(ax)
    improvePlot
    f3.Position = [100 100 1410 600];

    % comparing all modules
    
    t_struct = time_multiplier_singlebody * extract_runtime(profile_singlebody,'structures_one_case');
    t_geom   = time_multiplier_singlebody * extract_runtime(profile_singlebody,'geometry');
    t_econ   = time_multiplier_singlebody * extract_runtime(profile_singlebody,'econ');

    t_modules = [t_geom,t_hydro,t_multibody,t_struct,t_econ];
    name_modules = {'Geometry','Hydrodynamics','Dynamics (2-DOF)','Structures','Economics'};
    name_modules = reordercats(categorical(name_modules),name_modules);

    f4 = figure;
    h = bar(name_modules,t_modules);
    if ~isMATLABReleaseOlderThan('R2024b')
        h(1).Labels = h(1).YData;
    end
    ylabel('Runtime (s)')
    title('Module Time Breakdown')
    split_labels_by_space(h.Parent)
    improvePlot
end

function runtime = extract_runtime(profile_results, func_name, line_num)

    % for inner funcs (with '>'), need to just take after the '>'
    func_list = {profile_results.FunctionTable.FunctionName};
    func_list(contains(func_list,'>')) = extractAfter( func_list(contains(func_list,'>')), '>');

    idx_func = strcmp(func_list, func_name);
    if ~exist('line_num','var') % all lines
        runtime = sum([profile_results.FunctionTable(idx_func).TotalTime]);
    else % specific lines
        line_matrix = profile_results.FunctionTable(idx_func).ExecutedLines;
        idx_lines = ismember(line_matrix(:,1),line_num);
        runtime = sum(line_matrix(idx_lines,3));
    end

end

function lineNumbers = findStringInFile(fileName, searchString)
% generated by matlab playground AI
    % Open the file for reading
    fid = fopen([fileName '.m'], 'r');
    if fid == -1
        error('File cannot be opened: %s', fileName);
    end

    lineNumbers = []; % Initialize an empty array to hold line numbers
    lineIndex = 1;    % Initialize line index

    % Read the file line by line
    while ~feof(fid)
        line = fgetl(fid); % Get the current line
        if contains(line, searchString) % Check if the line contains the search string
            lineNumbers(end + 1) = lineIndex; % Store the line number
        end
        lineIndex = lineIndex + 1; % Increment line index
    end

    % Close the file
    fclose(fid);
end

function split_labels_by_space(ax)
    old_labels = ax.YTickLabel;
    num_rows = 2;
    labelArray = cell([num_rows,length(old_labels)]);
    for i=1:length(old_labels)
        tmp = split(ax.YTickLabel(i));
        buffer = repmat({''}, num_rows - length(tmp), 1);
        labelArray(:,i) = [tmp buffer];
    end
    pattern = ['%s' repmat('\\newline%s' ,[1 num_rows-1]) '\n'];
    new_labels = strtrim(sprintf(pattern, labelArray{:}));
    ax.YTickLabel = new_labels;
end
