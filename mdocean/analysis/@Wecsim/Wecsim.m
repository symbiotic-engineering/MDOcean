classdef Wecsim < GenericAnalysis
    %WECSIM Analysis class for WEC-Sim validation figures and tables
    %   Generates WEC-Sim error histograms/contour plots and validation tables

    properties
        fig_names = fig_names_wrapper();
        
        tab_names = {'WECSim_errors'};
    end

    methods (Static)
        intermed_result_struct = analysis_fcn(p,b)

        [fig_array,...
         tab_array_display,...
         tab_array_latex,...
         end_result_struct,...
         tab_firstrows,...
         tab_colspecs] = post_process_fcn(intermed_result_struct)
    end
end

function fig_names = fig_names_wrapper()
    [~,case_desc_cell] = make_all_cases(false);

    % fixme this is defined separately in analysis_fcn
    histogram_sep_figures = false;

    % fixme this is defined separately in power_matrix_compare
    contour_names = {'power_mech_unsat','power_elec_sat','CW_to_CW_max',...
        'float_amplitude','relative_amplitude','spar_amplitude','PTO_damping',...
        'force_pto','float_drag_force_fund','spar_drag_force_fund',...
        'float_drag_force_phase','spar_drag_force_phase',...
        'float_phase','spar_phase','rel_phase'};

    if histogram_sep_figures
        plot_names = ['histogram',contour_names];
        case_names = extractAfter([case_desc_cell{:}],'wecsim_');
        [case_names_mesh, plot_names_mesh] = meshgrid(case_names, plot_names);
        fig_names = strcat('wecsim_', case_names_mesh(:), '__', plot_names_mesh(:));
    
    else
        group_names = cell(size(case_desc_cell)); % 1 x 3
        num_groups = length(case_desc_cell); % 3
        num_cases_per_group = length(case_desc_cell{1}); % 4
        num_contours_per_case = length(contour_names); % 15
        plot_names = cell(num_contours_per_case * num_cases_per_group + 1,  num_groups); % (4*15+1=61) x 3

        for i = 1:num_groups % create separate plot names for each group
            case_names_i = extractAfter(case_desc_cell{i},'wecsim_');
            group_names{i} = extractBefore(case_names_i{1},'_drag');
            [case_names_mesh_i, contour_names_mesh_i] = meshgrid(case_names_i, contour_names);
            contour_plot_names_i = strcat(case_names_mesh_i(:), '__', contour_names_mesh_i(:)); % 60 per group
            plot_names(1:end-1,i) = contour_plot_names_i;
        end
        histogram_names = strcat(group_names,'__histogram');
        plot_names(end,:) = histogram_names;
        fig_names = strcat('wecsim_',plot_names(:)).';
    end
end