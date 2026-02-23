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
results = NaN(n,4);
timeout_sec = 10;

%% Loop through commits
for k = 1:10:61
    hash = hash_cell{k};
    fprintf('Checking out commit %d/%d: %s\n', k, n, hash);

    % Checkout commit
    system(sprintf('git checkout %s', hash));

    clear functions
    rehash

    try

        % Run asynchronously on worker
        f = parfeval(@()hydro_coeff_err(false), 1);
    
        completed = wait(f, 'finished', timeout_sec);
    
        if completed
            val = fetchOutputs(f);
        else
            fprintf('Timeout at commit %s\n', hash);
            cancel(f);
            val = NaN(1,4);
        end

        if numel(val) == 3
            val = [val NaN];   % pad with NaN
        end

        results(k,:) = val;

    catch ME
        warning('Error at commit %s: %s', hash, ME.message);
        results(k,:) = NaN;
    end
end

%% Restore original branch
fprintf('Restoring branch: %s\n', current_branch);
system(sprintf('git checkout %s', current_branch));

%% plot
names = {'A','B','|\gamma|','\angle \gamma'};
figure
for i=1:4
    if i==1 || i==3
        yyaxis left
    else
        yyaxis right
    end
    if i==1 || i==2
        linespec = '*-';
    else
        linespec = 'o';
    end
    plot(1:n, results(:,i), linespec,'DisplayName',names{i})
    hold on
end
xlabel(['Commits after ' start_date])
legend
