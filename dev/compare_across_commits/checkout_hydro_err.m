%% User settings
start_date = '2024-05-23';
end_date   = '2025-01-15';

%% Get current branch so we can restore it later
[~, current_branch] = system('git rev-parse --abbrev-ref HEAD');
current_branch = strtrim(current_branch);

%% Get commit hashes in date range (oldest → newest)
git_cmd = sprintf('git rev-list --reverse --since="%s" --until="%s" HEAD', ...
    start_date, end_date);

[~, commit_list_str] = system(git_cmd);
commit_list_cell = splitlines(commit_list_str);
commit_list_cell = commit_list_cell(2:2:end-1);

hash_cell = cellfun(@(s) s(1:min(7, end)), commit_list_cell, 'UniformOutput', false);
date_cell = cellfun(@(s) s(8:end),         commit_list_cell, 'UniformOutput', false);

n = numel(hash_cell);

%% Preallocate results
results = cell(n,1);
commit_ids = strings(n,1);

%% Loop through commits
for k = 1:n
    hash = hash_cell{k};
    fprintf('Checking out commit %d/%d: %s\n', k, n, hash);

    % Checkout commit
    system(sprintf('git checkout %s', hash));

    try
        % Clear function cache in case file changed
        clear hydro_coeff_err

        % Run your function
        narg = nargin('hydro_coeff_err');
        if narg==1
            val = hydro_coeff_err(false);
        elseif narg==0
            val = hydro_coeff_err;
        end

        results{k} = val;
        commit_ids(k) = hash;

    catch ME
        warning('Error at commit %s: %s', hash, ME.message);
        results{k} = NaN;
        commit_ids(k) = hash;
    end
end

%% Restore original branch
fprintf('Restoring branch: %s\n', current_branch);
system(sprintf('git checkout %s', current_branch));


%% Convert results to matrix if consistent size
try
    results_mat = cell2mat(results);
catch
    warning('Results have inconsistent sizes; kept as cell array.');
    results_mat = results;
end
