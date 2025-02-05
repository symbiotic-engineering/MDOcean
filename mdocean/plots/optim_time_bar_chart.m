
function fig = optim_time_bar_chart(suite,results)

full_names = {results.Name};
times = [results.Duration];

idx_tabs = contains(full_names,'which_figs=none');
idx_figs = contains(full_names,'which_tabs=none');

idx_matlab_fig = false(size(idx_figs));
for i=find(idx_figs)
    try
        fig_path = results(i).Details.DiagnosticRecord.TestDiagnosticResults.Artifacts(1).FullPath;
        fig = openfig(fig_path,'invisible');
        idx_matlab_fig(i) = isempty(fig.UserData);
    catch
        warning("Figure %s errored so it isn't included in the time comparison.\n",full_names{i})
    end
    
end

idx_matlab = idx_tabs | idx_matlab_fig;

fig_tab_names = cell(size(suite));

% table names - use even indices to fetch which_tab
st = suite(idx_tabs);
stp = [st.Parameterization];
if ~isempty(stp)
    fig_tab_names(idx_tabs) = {stp(2:2:length(stp)).Name};
end

% fig names - use odd indices to fetch which_fig
sf = suite(idx_matlab_fig);
sfp = [sf.Parameterization];
if ~isempty(sfp)
    fig_tab_names(idx_matlab_fig) = {sfp(1:2:length(sfp)).Name};
end

figure
fig = bar(categorical(remove_underscores(fig_tab_names(idx_matlab))),times(idx_matlab));
title('Runtime (seconds)')
improvePlot

end