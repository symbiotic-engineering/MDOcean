function intermed_result_struct = analysis_fcn(~,b)
    % Get constraints information
    % Store results for post-processing
    tab = cell2table(remove_underscores(b.constraint_names.'));
    rows = 1:24; % Only show slamming constraint for one sea state
    intermed_result_struct.constraint_names = tab(rows,:);
end
