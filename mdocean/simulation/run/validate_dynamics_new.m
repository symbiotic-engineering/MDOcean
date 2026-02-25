% % create cases
% 1 - singlebody, wamit geometry
% 1a - drag off, wamit coeffs
% 1b - drag on, wamit coeffs
% 1c - drag off, meem coeffs
% 1c2 - drag off, meem coeffs, N=50
% 1d - drag on, meem coeffs

% 2 - multibody, wamit geometry
% 2a - drag off, wamit coeffs
% 2b - drag on, wamit coeffs
% 2c - drag off, meem coeffs
% 2c2 - drag off, meem coeffs, N=50
% 2d - drag on, meem coeffs

% 3 - report geometry
% 3a - drag on, meem coeffs
% 3b - drag on, wamit coeffs

% % run cases
% all use run_dynamic_validation(X,p,reportOn)

% % plot cases
% all get power_matrix_compare + its report
% wecsim also gets a histogram

% % plot comparison between cases
% the histograms are put in the same figure
% and overall errors are put in a struct

% [errors_singlebody, errors_multibody, ...
%  errors_report, T, fig_singlebody, fig_multibody] = validate_dynamics()

function [case_cell, filename_cell, err_structs, fig_mats] = validate()
    % override to have fewer sea states for the sake of fast debugging
    runOnlyFewSeaStates = false;

    case_cell = make_all_cases(runOnlyFewSeaStates);
    filename_cell = run_wecsim_all_cases(case_cell);
    err_struct_figs_cell = plot_all_cases(case_cell,filename_cell,runOnlyFewSeaStates);
    err_structs = err_struct_figs_cell{:,1};
    fig_mats = err_struct_figs_cell{:,2};

    errors_wecsim_sb = err_structs{1};
    errors_wecsim_mb = err_structs{2};
    errors_report_mb = err_structs{3}; % first index is wecsim to report error, second index is mdocean to report error

    %% make table comparing error of mdocean to various ground truths
    errors_report_mdocean = structfun(@(s)s(2), errors_report_mb,'UniformOutput',false); % just mdocean to report error
    T_report = struct2table(errors_report_mdocean);
    T_singlebody = struct2table(errors_wecsim_sb);
    T_multibody = struct2table(errors_wecsim_mb);
    T = [T_report;T_singlebody;T_multibody];
    T.Row = {'RM3 Report','WEC-Sim Singlebody','WEC-Sim Multibody'};
end

function out_cell = run_wecsim_all_cases(case_cell)
    num_case_groups = size(case_cell,1);
    out_cell = cell(1,num_case_groups);
    for case_group_idx = 1:num_case_groups
        X = case_cell{case_group_idx,1};
        p_array = case_cell{case_group_idx,2};
        out_each_case_array = cell(1,numel(p_array));
        for p_idx = 1:numel(p_array)
            p = p_array(p_idx);
            out_each_case_array(p_idx) = run_wecsim_validation(X,p);
        end
        out_cell(case_group_idx) = out_each_case_array;
    end
end

function output_filename = run_wecsim_validation(X,p)
    % X and p need to be in the workspace for runRM3Parallel script to work right
    runRM3Parallel % this script uses p and modifies it, and saves output_filename to workspace
end

function err_struct_figs_cell = plot_all_cases(case_cell,filename_cell,runOnlyFewSeaStates)
    histogram_sep_figures = true;
    num_case_groups = size(case_cell,1);
    err_struct_figs_cell = cell(num_case_groups,2);

    for case_group_idx = 1:num_case_groups
        X = case_cell{case_group_idx,1};
        p_array = case_cell{case_group_idx,2};
        filename_array = filename_cell{case_group_idx};
        out_each_case_array = cell(1,numel(p_array));

        if ~histogram_sep_figures
            hist_fig = figure;
            t = tiledlayout(1, length(p_array));
        end


        my_fieldnames = ["pct_error_baseline",...
                        "pct_error_drag",...
                        "pct_error_total"];
        my_subtitles = ["Drag Off, Identical Hydro Coefficients",...
                        "Drag On, Identical Hydro Coefficients",...
                        "Drag On, Different Hydro Coefficients"];

        for p_idx = 1:numel(p_array)
            p = p_array(p_idx);
            filename = filename_array(p_idx);
            if p.use_multibody
                widths = [.15 2 5];
            else
                widths = [.05 .25 5];
            end
            width = widths(p_idx);

            if histogram_sep_figures
                hist_fig = figure;
                ax = gca;
            else
                figure(hist_fig);
                ax = nexttile(t);
                subtitle(my_subtitles(p_idx),'FontSize',13)
            end
            [weighted_pwr_err,...
            max_amp_err,...
            pwr_err,amp_err,fig_vec] = plot_per_case(X,p,filename,RM3reportOn,...
                                                    runOnlyFewSeaStates,ax,width);

            if histogram_sep_figures
                fig_vec = [fig_vec hist_fig];
            end
            
            error_struct.(my_fieldnames(p_idx)) = [weighted_pwr_err, max_amp_err];
            % weighted_pwr_err_array(p_idx) = weighted_pwr_err;
            % max_amp_err_array(p_idx) = max_amp_err;
            % pwr_err_array(p_idx) = pwr_err;
            % amp_err_array(p_idx) = amp_err;
            fig_mat(p_idx,:) = fig_vec;
        end
        if ~histogram_sep_figures
            leg = legend;
            title(leg,'Error in:')
            leg.Position = [0.514,0.212,0.226,0.222];
            
            if multibody
                mb_string = 'Multi-Body';
            else
                mb_string = 'Single Body';
            end
            xlabel(t, 'Percent Error of MDOcean to WecSim')
            ylabel(t, 'Fraction of Sea States')
            title(t,[mb_string ' Dynamics'])
            improvePlot
            hist_fig.Position = [7.4,137,1523.2,600];
        end
        err_struct_figs_cell(case_group_idx,:) = {error_struct, fig_mat};
    end
end

function [weighted_pwr_err,...
          max_amp_err,...
          pwr_err,amp_err,figs] = plot_per_case(X,p,wecsim_filename,RM3reportOn,...
                                                runOnlyFewSeaStates,ax,width)

    [weighted_pwr_err, max_amp_err, ...
     pwr_err,          amp_err,    figs] = power_matrix_compare(X,p,wecsim_filename, ...
                                                    RM3reportOn,runOnlyFewSeaStates);

    make_report(figs,wecsim_filename,p)

    make_histogram_on_axis(ax,width,...
                            pwr_err,amp_err,...
                            weighted_pwr_err,max_amp_err)

end

function make_report(figs,wecsim_filename,p)
    import mlreportgen.report.* 
    import mlreportgen.dom.* 
    rpt = Report(['DynamicValidation_' wecsim_filename],'pdf'); 
    rpt.Layout.Landscape = true;
    pm = PageMargins();
    pm.Left='0.5in'; pm.Right='0.5in';
    rpt.Layout.PageMargins = pm;
    
    % add parameters
    display_fields = {'use_MEEM','use_multibody','C_d_float','C_d_spar','harmonics'};
    for i=1:length(display_fields)
        pDisplay.(display_fields{i}) = p.(display_fields{i});
    end
    append(rpt,struct2table(pDisplay))

    % add figures
    for i=1:length(figs)
        fig = figs(i);
        pos = fig.Position;
        fig.Position = [pos(1) pos(2) 2*pos(3) pos(4)];
        if isvalid(fig)
            try
                append(rpt,Figure(fig))
            catch
                warning(['Figure %d of %d has been deleted and was skipped ' ...
                'from the report.'],i,length(figs))
            end
        else
            warning(['Figure %d of %d has been deleted and was skipped ' ...
                'from the report.'],i,length(figs))
        end
    end

    close(rpt)
    %rptview(rpt)
end

function make_histogram_on_axis(ax,width,...
                                pwr_err,amp_err,...
                                weighted_pwr_err,max_amp_err)
    % box plot
    axes(ax)
    histogram(ax,pwr_err(:),'Normalization','probability','BinWidth',width,'DisplayName','Mechanical Power')
    hold(ax,'on')
    h = histogram(ax,amp_err(:),'Normalization','probability','BinWidth',width,'HandleVisibility','off');
    hb = bar(ax,h.BinEdges(1:end-1),-h.Values,'histc');
    hb.FaceColor = [0.8500 0.3250 0.0980];
    hb.FaceAlpha = 0.6; 
    hb.DisplayName = 'Float Amplitude';

    h.EdgeAlpha = 0; h.FaceAlpha = 0; % transparent

    ylim([-.55 .55])
    xx = xlim;
    xlim([-1 1]*max(abs(xx))) % zero centered on x
    yy = ylim;

    plot(ax,[1 1]*weighted_pwr_err,[0 yy(2)],'Color',[0 0.4470 0.7410], ...
        'DisplayName','JPD-Weighted Average Power')
    plot(ax,[1 1]*max_amp_err,     [yy(1) 0],'Color',[0.8500 0.3250 0.0980], ...
        'DisplayName','Maximum Float Amplitude')

    text_offset_if_neg = -max(abs(xx))/3.7;
    text_offset_if_pos =  max(abs(xx))/18;
    x_text_pwr = weighted_pwr_err + text_offset_if_neg*logical(weighted_pwr_err<0) ...
                + text_offset_if_pos*logical(weighted_pwr_err>0);
    x_text_amp = max_amp_err      + text_offset_if_neg*logical(max_amp_err<0) ...
                + text_offset_if_pos*logical(max_amp_err>0);

    text(x_text_pwr, yy(2)*.9, sprintf('%+0.1f%%',weighted_pwr_err), ...
        "Color",[0 0.4470 0.7410],'FontSize',12)
    text(x_text_amp, yy(1)*.9, sprintf('%+0.1f%%',max_amp_err), ...
        "Color",[0.8500 0.3250 0.0980],'FontSize',12)

    plot(ax,[0 0],                  yy, 'k--','HandleVisibility','off')
    plot(ax,[-1 1]*max(abs(xx)),[0 0],'k-','LineWidth',.5,'HandleVisibility','off')
    xtickformat('percentage')

end
