function [errors_singlebody, errors_multibody, ...
            errors_report, T, fig_singlebody, fig_multibody] = validate_dynamics()
    %% Round 1: singlebody, wamit geometry
    [errors_singlebody,fig_singlebody] = wecsim_error_breakdown(false);
    
    %% Round 2: multibody, wamit geometry
    [errors_multibody,fig_multibody] = wecsim_error_breakdown(true);
    
    %% Round 3: report geometry
    errors_report = report_error_breakdown(); % first index is wecsim to report error, second index is mdocean to report error

    %% make table comparing error of mdocean to various ground truths
    errors_report_mdocean = structfun(@(s)s(2), errors_report,'UniformOutput',false); % just mdocean to report error
    T_report = struct2table(errors_report_mdocean);
    T_singlebody = struct2table(errors_singlebody);
    T_multibody = struct2table(errors_multibody);
    T = [T_report;T_singlebody;T_multibody];
    T.Row = {'RM3 Report','WEC-Sim Singlebody','WEC-Sim Multibody'};

end

% bonus that I have't checked yet: sweep drag, irregular waves, force saturation


function [weighted_pwr_err,max_amp_err,pwr_err,amp_err] = run_dynamic_validation(X,p,RM3reportOn)
    if nargin<3
        RM3reportOn = false;
    end

    % override to have fewer sea states for the sake of fast debugging
    runOnlyFewSeaStates = true;
    if runOnlyFewSeaStates
        p.Hs = p.Hs(3:4);
        p.T = p.T(4:5);
        p.JPD = p.JPD(3:4,4:5);
    end
    % rerun wecsim or load files
    runWecSim = true;
    if ~runWecSim
        if ~p.use_MEEM
            meem = 'off';
        else
            meem = num2str(p.harmonics);
        end
        C_d_s = num2str(p.C_d_spar);
        C_d_f = num2str(p.C_d_float);
        mb = num2str(p.use_multibody);
        
        wecsim_filename_start = ['wecsim_sparcd' C_d_s '_floatcd' C_d_f '_multibody_' mb ...
                            '_meem_' meem];
        d = dir(['results_3_31\**\' wecsim_filename_start '*'] );
        wecsim_filename = d.name;
    else
        wecsim_filename = run_wecsim_validation(X,p);
    end
    [weighted_pwr_err,max_amp_err,pwr_err,amp_err,figs] = power_matrix_compare(X,p,wecsim_filename, ...
                                                    RM3reportOn,runOnlyFewSeaStates);

    make_report(figs,wecsim_filename,p)
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

function output_filename = run_wecsim_validation(X,p)
    % X and p need to be in the workspace for runRM3Parallel script to work right
    runRM3Parallel % this script uses p and modifies it, and saves output_filename to workspace

end

function [errors,f] = wecsim_error_breakdown(multibody)
    f = figure;
    t = tiledlayout(1, 3);

    b = var_bounds('wecsim');
    X = [b.X_noms; 1];
    X(7) = 1e8; % disable power limit
    
    % 1. drag off, wamit coeffs, wecsim geometry: should match wecsim very well <2%
    p = parameters('wecsim');
    p.use_multibody = multibody;
    p.C_d_float = 0;
    p.C_d_spar = 0;
    p.use_MEEM = false;
    X(strcmp(b.var_names,'F_max')) = Inf;
    ax1 = nexttile(t);
    errors.pct_error_baseline = create_wecsim_validation_histogram(X,p,ax1,.15);
    subtitle('Drag Off, Identical Hydro Coefficients','FontSize',13)
    
    % 2. 1 but drag on: gives me % error that comes from drag
    p = parameters('wecsim');
    p.use_multibody = multibody;
    p.use_MEEM = false;
    X(strcmp(b.var_names,'F_max')) = Inf;
    figure(f);
    ax2 = nexttile(t);
    errors.pct_error_drag = create_wecsim_validation_histogram(X,p,ax2,.25);
    subtitle('Drag On, Identical Hydro Coefficients','FontSize',13)
%     
%     % 3. drag back off but meem coffs: gives % error that comes from meem
%     p = parameters('wecsim');
%     p.use_multibody = multibody;
%     p.C_d_float = 0;
%     p.C_d_spar = 0;
%     errors.pct_error_meem_all = run_dynamic_validation(X,p);
%     
%     %     3a. meem with N=50: gives meem truncation error
%     p = parameters('wecsim');
%     p.use_multibody = multibody;
%     p.C_d_float = 0;
%     p.C_d_spar = 0;
%     p.harmonics = 50;
%     errors.pct_error_meem_all = run_dynamic_validation(X,p);
    
    %     3b. wamit with zero phase: gives meem phase error
    % right now this requires manual fiddling with wamit coeffs
    
    % 4. drag meem interaction
    p = parameters('wecsim');
    p.use_multibody = multibody;
    X(strcmp(b.var_names,'F_max')) = Inf;
    figure(f);
    ax3 = nexttile(3);
    errors.pct_error_total = create_wecsim_validation_histogram(X,p,ax3,5);
    subtitle('Drag On, Different Hydro Coefficients','FontSize',13)
    leg = legend;
    title(leg,'Error in:')
    leg.Position = [0.537,0.212,0.226,0.222];
    
    if multibody
        mb_string = 'Multi-Body';
    else
        mb_string = 'Single Body';
    end
    xlabel(t, 'Percent Error of MDOcean to WecSim')
    ylabel(t, 'Fraction of Sea States')
    title(t,[mb_string ' Dynamics'])
    improvePlot
    f.Position = [7.4,137,1523.2,600];
end

function [weighted_pwr_err,max_amp_err] = create_wecsim_validation_histogram(X,p,ax,width)
    [weighted_pwr_err,...
     max_amp_err,...
     pwr_err,amp_err] = run_dynamic_validation(X,p);
    
    if nargin>2
        % box plot
        axes(ax)
        histogram(ax,pwr_err(:),'Normalization','probability','BinWidth',width,'DisplayName','Mechanical Power')
        hold(ax,'on')
        h = histogram(ax,amp_err(:),'Normalization','probability','BinWidth',width,'HandleVisibility','off');
        hb = bar(ax,h.BinEdges(1:end-1),-h.Values,'histc');
        hb.FaceColor = [0.8500 0.3250 0.0980];
        hb.FaceAlpha = h.FaceAlpha; 
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
end

function errors = report_error_breakdown()
    p = parameters();
    b = var_bounds();
    X = [b.X_noms; 1];
    
    % 1. drag on, meem coeffs, report geometry, try to match report 10%?
    errors.pct_error_total = run_dynamic_validation(X,p,true);
    
    % 2. 1 but wamit coeffs and multibody
    p.use_multibody = true;
    p.use_MEEM = false;
    X(strcmp(b.var_names,'F_max')) = Inf;
    errors.pct_error_baseline = run_dynamic_validation(X,p,true);

    % drag: placeholder for now
    errors.pct_error_drag = [NaN NaN];
end