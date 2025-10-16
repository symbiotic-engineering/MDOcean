function move_results_to_pubs()

    [figs_in_RE, figs_in_AOR, tabs_in_RE, tabs_in_AOR] = fig_tab_pub_mapping();

    RE = 'pubs/renewable-energy-mdo';
    AOR = 'pubs/applied-ocean-research-model';

    success = true;

    for i = 1:length(figs_in_RE)
        success = success && copy_result(figs_in_RE{i}, '.pdf', RE);
    end
    for i = 1:length(figs_in_AOR)
        success = success && copy_result(figs_in_AOR{i}, '.pdf', AOR);
    end
    for i = 1:length(tabs_in_RE)
        success = success && copy_result(tabs_in_RE{i}, '.tex', RE);
    end
    for i = 1:length(tabs_in_AOR)
        success = success && copy_result(tabs_in_AOR{i}, '.tex', AOR);
    end
    assert(success, 'Some files were not found and could not be copied.')
end

function success = copy_result(class_dot_description, extension, dest_folder)
    s = split(class_dot_description, '.');
    class_name = s{1};
    fig_tab_desc = s{2};
    src = fullfile('results', class_name, [fig_tab_desc, extension]);
    dest = fullfile(dest_folder, [fig_tab_desc, extension]);
    if exist(src,"file")
        copyfile(src, dest);
        success = true;
    else
        warning(['File not found: ' src])
        success = false;
    end
end