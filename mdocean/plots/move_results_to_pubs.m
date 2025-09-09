function move_results_to_pubs()

    [figs_in_RE, figs_in_AOR, tabs_in_RE, tabs_in_AOR] = fig_tab_pub_mapping();

    RE = 'pubs/renewable-energy-mdo';
    AOR = 'pubs/applied-ocean-research-model';

    for i = 1:length(figs_in_RE)
        copy_result(figs_in_RE{i}, '.pdf', RE);
    end
    for i = 1:length(figs_in_AOR)
        copy_result(figs_in_AOR{i}, '.pdf', AOR);
    end
    for i = 1:length(tabs_in_RE)
        copy_result(tabs_in_RE{i}, '.tex', RE);
    end
    for i = 1:length(tabs_in_AOR)
        copy_result(tabs_in_AOR{i}, '.tex', AOR);
    end

end

function copy_result(class_dot_description, extension, dest_folder)
    s = split(class_dot_description, '.');
    class_name = s{1};
    fig_tab_desc = s{2};
    src = fullfile('results', class_name, [fig_tab_desc, extension]);
    dest = fullfile(dest_folder, [fig_tab_desc, extension]);
    copyfile(src, dest);
end