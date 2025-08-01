function figOut = unittest_time_bar_chart(suite,results)

full_names = {results.Name};
times = [results.Duration];

idx_tabs = contains(full_names,'which_figs=none');
idx_figs = contains(full_names,'which_tabs=none');

idx_matlab_fig = false(size(idx_figs));
for i=find(idx_figs)
    try
        tdr = [results(i).Details.DiagnosticRecord.TestDiagnosticResults];
        fig_path = tdr.Artifacts(1).FullPath;
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

figOut = figure;
hold on
bar(categorical(remove_underscores(fig_tab_names(idx_matlab))),times(idx_matlab));
for i = 1:length(idx_matlab)
    if times(i) == 0
        plot(i,0,'x')
    end
end
title('Runtime (seconds)')
hold off
improvePlot

end