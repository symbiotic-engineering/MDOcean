function [figs, runtimeLocal, runtimeGlobal] = param_sweep(filename_uuid)

    %% Setup
    b = var_bounds();
    if nargin>0
        b.filename_uuid = filename_uuid;
    end
    [p,T] = parameters();
    
    param_names = T.name_pretty(T.sweep);  % list of parameters to sweep
    params = T.name(T.sweep);
    p_val  = [T.value_normalize{T.sweep}];
    p_idxs = T.index_normalize(T.sweep);

    dvar_names = b.var_names_pretty(1:end-1);
    
    groups = categorical(T.subsystem(T.sweep));
    color_groupings = {'b','g','r','k','y','m'};
    colors = color_groupings(groups);
    
    %%
    % use the optimal x as x0 to speed up the sweeps
    % and obtain gradients
    x0 = b.X_start_struct;
    x0_vec_1 = gradient_optim(x0,p,b);
    x0_struct_1 = cell2struct(num2cell(x0_vec_1(1:end-1,:)),b.var_names(1:end-1)',1);

    % rerun nominal optim to check that it's the same as the first time
    [x0_vec, J0, ~, ~, lambdas, grads, hesses] = gradient_optim(x0_struct_1,p,b);
    same_as_before = ismembertol(x0_vec,x0_vec_1,1e-6);
    if ~all(same_as_before,'all')
        str = append(['First optimization of nominal parameters gave different x* than second. Using second.' newline ...
                      '#1 x*:'], newline, formattedDisplayText(x0_vec_1.'),'#2 x*:' ,newline, formattedDisplayText(x0_vec.'), 'Δ x*:', newline, formattedDisplayText(x0_vec_1.'-x0_vec.'));
        warning(str)
    end

    num_constr_nl = length(b.constraint_names);
    g_lambda_0(:,1) = combine_g_lambda(lambdas(1),x0_vec(:,1),p,b);
    g_lambda_0(:,2) = combine_g_lambda(lambdas(2),x0_vec(:,2),p,b);
    
    %% Obtain normalized local sensitivity 
    disp('done optimizing, starting local param sensitivity')
    t = tic;
    [par_x_star_par_p_local, ...
     dJ_star_dp_lin_local, dJ_star_dp_quad_local, ...
     par_J_par_p_local, ...
     delta_p_change_activity_local] = local_sens_both_obj_all_param(x0_vec, J0, p, params, p_val, p_idxs, ...
                                                                    lambdas, grads, hesses, num_constr_nl);
    runtimeLocal = toc(t);

    %% Obtain normalized global sensitivity
    disp('done local, starting global param sensitivity')
    t = tic;
    [par_x_star_par_p_global, dJstar_dp_global, ...
     delta_p_change_activity_global,...
     ~, ~, ~, ~, figs_global] = global_sens_all_param(params, param_names, p_val, ...
                                                                x0_vec, dvar_names, J0, g_lambda_0, ...
                                                                p, b, colors, groups);
    runtimeGlobal = toc(t);

    %% Post processing
    dJdp_combined = [par_J_par_p_local; dJ_star_dp_quad_local; dJ_star_dp_lin_local; dJstar_dp_global];
    dJdp_names = {'partial','total linear','total quadratic','re-optimization'};
    
    % color grid plots
    % fig 1: dJ*/dp combined
    f1 = sensitivity_plot(dJdp_combined, 'dJ*/dp normalized', param_names, dJdp_names, ...
        'Parameters p', 'Type of Sensitivity');
    % fig 2: dJ*/dp global
    f2 = sensitivity_plot(dJstar_dp_global, 'dJ*/dp normalized: global', param_names, '', ...
        'Parameters p', '');
    
    % fig 3: dx*/dp local
    f3 = dxdp_plot(par_x_star_par_p_local,  param_names, dvar_names, 'local');
    % fig 4: dx*/dp global
    f4 = dxdp_plot(par_x_star_par_p_global, param_names, dvar_names, 'global');
    
    % fig 5: delta p local
    f5 = delta_p_plot(delta_p_change_activity_local,  b, dvar_names, param_names, 'local');
    % fig 6: delta p global
    f6 = delta_p_plot(delta_p_change_activity_global, b, dvar_names, param_names, 'global');

    figs = [figs_global,f1,f2,f3,f4,f5,f6]; % figs_global has 6+num_dvs=18 so 24 total

end

function fig = delta_p_plot(delta_p_norm,b,dvar_names,param_names,title_suffix)
    % combine all the slamming constraints into one so it fits on plot
    constr_names  = [b.constraint_names_pretty b.lin_constraint_names_pretty ...
        strcat(dvar_names," lower"), strcat(dvar_names," upper")];
    idx_slam = contains(constr_names,'Slamming');
    delta_p_norm_combine_slamming = delta_p_norm(:,~idx_slam);
    delta_p_norm_slam = delta_p_norm(:,idx_slam);
    [~,idx_delta_p_slam] = min(abs(delta_p_norm_slam),[],2);
    idx_delta_p_slam = sub2ind(size(delta_p_norm_slam),1:length(param_names),idx_delta_p_slam');
    delta_p_norm_combine_slamming(:,end+1) = delta_p_norm_slam(idx_delta_p_slam);
    
    titl = ['Normalized \Deltap to change constraint activity: ' title_suffix];
    xticks = [constr_names(~idx_slam),'Slamming'];
    yticks = param_names;
    xlab = 'Constraints';
    ylab = 'Parameters';
    fig = sensitivity_plot(delta_p_norm_combine_slamming, titl, xticks, yticks, xlab, ylab);
end

function fig = dxdp_plot(dxdp, param_names, dvar_names,title_suffix)
    titl = ['dx*/dp normalized: ' title_suffix];
    xticks = param_names;
    yticks = dvar_names;
    xlab = 'Parameters p';
    ylab = 'Design Variables x';
    fig = sensitivity_plot(dxdp, titl, xticks, yticks, xlab, ylab);
end

function fig = sensitivity_plot(matrix, titl, xticks, yticks, xlab, ylab)
    fig = color_each_element(matrix);
    ax = gca;
    title(titl)
    set(ax,'XTickLabel',xticks)
    set(ax,'YTickLabel',yticks)
    clim([-3 3])
    colormap(bluewhitered)
    xlabel(xlab)
    ylabel(ylab)
    improvePlot
    set(ax,'XTickLabelRotation', 60)
    if length(xticks) < 50
        ax.XAxis.FontSize = 14;
    else
        ax.XAxis.FontSize = 10;
    end
    if length(yticks) > 40
        ax.YAxis.FontSize = 10;
    end
    set(fig,'Position',[1 41 1536 840]) % full screen
end

%% Rerun optimization for global sensitivities
function [par_x_star_par_p_global, ...
          dJstar_dp_global, ...
          delta_p_change_activity_global,...
          LCOE, P_var, X_LCOE, X_Pvar, figs] = global_sens_all_param(params, param_names, p_val, ...
                                                               x0_vec, dvar_names, J0, g_lambda_0, ...
                                                               p, b, colors, groups)
    %ratios = .8 : .1 : 1.2;
    ratios = [.9 .98 1 1.02 1.1];
    num_DVs = size(x0_vec,1)-1;
    num_constr = length(b.constraint_names) + length(b.lin_constraint_names) + 2*num_DVs;
    
    [LCOE, P_var]    = deal(zeros(length(params), length(ratios)));
    [X_LCOE, X_Pvar] = deal(zeros(length(params), length(ratios), num_DVs));
    [g_lambda_LCOE, g_lambda_Pvar] = deal(zeros(length(params), length(ratios), num_constr));
    x0 = cell2struct(num2cell(x0_vec),b.var_names',1);

    parfor i=1:length(params)
        param_name = params{i};
        [LCOE(i,:),  X_LCOE(i,:,:), ...
         g_lambda_LCOE(i,:,:),...
         P_var(i,:), X_Pvar(i,:,:),...
         g_lambda_Pvar(i,:,:)] = global_sens_reoptimize(x0, p, b, param_name, ratios, num_DVs, num_constr);
    
    end
    
    % fill in ratios == 1 (left blank since same as initial)
    X_LCOE_0_3d = permute(repmat(x0_vec(1:end-1,1).',length(params),1,sum(ratios==1)),[1 3 2]);
    X_Pvar_0_3d = permute(repmat(x0_vec(1:end-1,2).',length(params),1,sum(ratios==1)),[1 3 2]);
    g_lambda_LCOE_0_3d = permute(repmat(g_lambda_0(:,1).',length(params),1,sum(ratios==1)),[1 3 2]);
    g_lambda_Pvar_0_3d = permute(repmat(g_lambda_0(:,2).',length(params),1,sum(ratios==1)),[1 3 2]);
    LCOE(:, ratios==1) = J0(1);
    P_var(:,ratios==1) = J0(2);
    X_LCOE(:,ratios==1,:) = X_LCOE_0_3d;
    X_Pvar(:,ratios==1,:) = X_Pvar_0_3d;
    g_lambda_LCOE(:,ratios==1,:) = g_lambda_LCOE_0_3d;
    g_lambda_Pvar(:,ratios==1,:) = g_lambda_Pvar_0_3d;
    
    % check if deltas too small, indicating potential finite precision error
    delta_x_LCOE = abs(X_LCOE - X_LCOE_0_3d(:,1,:));
    delta_x_Pvar = abs(X_Pvar - X_Pvar_0_3d(:,1,:));
    delta = [delta_x_LCOE delta_x_Pvar];
    delta_too_small = (delta < 1e-6) & (delta ~= 0);
    if any(delta_too_small,'all') 
        warning(['The global sensitivity is potentially inaccurate ' ...
            'due to finite precision effects. You must either decrease ' ...
            'options.StepTolerance in gradient_optim, or increase ' ...
            'abs(ratios-1) in param_sweep.'])
    end
    
    % info for optimal design with nominal parameter values (not nominal design)
    col_nom = find(ratios==1);
    LCOE_nom = LCOE(1,col_nom);
    Pvar_nom = P_var(1,col_nom);
    X_LCOE_nom = squeeze(X_LCOE(1,col_nom,:));
    X_Pvar_nom = squeeze(X_Pvar(1,col_nom,:));
    g_lambda_LCOE_nom = squeeze(g_lambda_LCOE(1,col_nom,:)).';
    g_lambda_Pvar_nom = squeeze(g_lambda_Pvar(1,col_nom,:)).';
    
    %% Calculate objective sensitivities
    
    slope_LCOE = get_slope(LCOE, ratios);
    slope_Pvar = get_slope(P_var, ratios);
    
    slope_LCOE_norm = slope_LCOE.' / LCOE_nom; % normalize
    slope_Pvar_norm = slope_Pvar.' / Pvar_nom;

    if all(~isfinite(slope_LCOE),'all') || all(~isfinite(slope_Pvar),'all')
        msg = ['All slopes are NaN, meaning all optimizations failed. ' ...
            'This might be because no feasible solution can be found. ' ...
            'slope_LCOE = ', num2str(slope_LCOE), ' and slope_Pvar = ', num2str(slope_Pvar)];
        error(msg)
    end
    
    %% Calculate design variable sensitivities
    slope_X_LCOE = get_slope(X_LCOE, ratios);
    slope_X_Pvar = get_slope(X_Pvar, ratios);
    
    slope_X_LCOE_norm  = slope_X_LCOE ./ X_LCOE_nom.'; % normalize
    slope_X_Pvar_norm  = slope_X_Pvar ./ X_Pvar_nom.';
    
    %% calculate lagrange multiplier and constraint sensitivities
    dlambda_g_dp_LCOE = get_slope(g_lambda_LCOE, ratios);
    dlambda_g_dp_Pvar = get_slope(g_lambda_Pvar, ratios);
    delta_p_LCOE = -g_lambda_LCOE_nom ./ dlambda_g_dp_LCOE;
    delta_p_Pvar = -g_lambda_Pvar_nom ./ dlambda_g_dp_Pvar;
    delta_p_LCOE_norm = delta_p_LCOE ./ p_val.';
    delta_p_Pvar_norm = delta_p_Pvar ./ p_val.';

    %% assign normalized outputs - fixme ignoring slope_Pvar for now
    dJstar_dp_global = slope_LCOE_norm;
    par_x_star_par_p_global = slope_X_LCOE_norm.';
    delta_p_change_activity_global = delta_p_LCOE_norm;
    
    %% Line plots showing nonlinearity
    fig_num_start = gcf().Number+1;
    f1 = figure(fig_num_start);
    subplot 121
    plot(ratios,LCOE/LCOE_nom)
    xlabel('Parameter Ratio from nominal')
    ylabel('LCOE ratio from nominal')
    legend(param_names)
    improvePlot
    grid on
    
    subplot 122
    plot(ratios,P_var/Pvar_nom)
    xlabel('Parameter ratio from nominal')
    ylabel('Power Variation ratio from nominal')
    legend(param_names)
    improvePlot
    grid on
    set(f1,'Position',[1 41 1536 840]) % full screen
    
    for i = 1:num_DVs
        f2 = figure(fig_num_start+1);
        subplot(3,5,i)
        plot(ratios, X_LCOE(:,:,i)./X_LCOE_nom(i))
        xlabel ('Parameter ratio from nominal')
        ylabel ('X* ratio from nominal')
        title([dvar_names{i} ' - min LCOE'])
        %improvePlot
        grid on
        
        f3 = figure(fig_num_start+2);
        subplot(3,5,i)
        plot(ratios, X_Pvar(:,:,i)./X_Pvar_nom(i))
        xlabel ('Parameter ratio from nominal')
        ylabel ('X* ratio from nominal')
        title([dvar_names{i} ' - min c_v'])
        %improvePlot
        grid on
    end
    legend(param_names)
    set(f2,'Position',[1 41 1536 840]) % full screen
    set(f3,'Position',[1 41 1536 840]) % full screen

    %% Tornado chart for overall slope
    f4 = sensitivity_tornado_barh(slope_LCOE_norm, slope_Pvar_norm, param_names, colors, groups, 'dJ*/dp');
    
    % separate figures for each design variable, with subplot for each objective
    for i=1:num_DVs
        f5_arr(i) = sensitivity_tornado_barh(slope_X_LCOE_norm(:,i), slope_X_Pvar_norm(:,i), ...
                                param_names, colors, groups, ['d' dvar_names{i} '^*/dp']);
    end

    % separate figures for each objective, with subplots for each design variable
    f6 = figure;
    t = tiledlayout(3,4);
    for i=1:num_DVs
        nexttile
        sensitivity_tornado_barh_inner(slope_X_LCOE_norm(:,i), param_names, colors, 8)
        title(['d' dvar_names{i} '^*/dp'])
    end
    color_legend(groups,colors,'eastoutside')
    title(t,'Normalized Sensitivities at Minimum LCOE')
    t.TileSpacing = 'compact';
    t.Padding = 'compact';

    f7 = figure;
    t2 = tiledlayout(3,4);
    for i=1:num_DVs
        nexttile
        sensitivity_tornado_barh_inner(slope_X_Pvar_norm(:,i), param_names, colors, 8)
        title(['d' dvar_names{i} '^*/dp'])
    end
    color_legend(groups,colors,'eastoutside')
    title(t2,'Normalized Sensitivities at Minimum Design Cost')
    t2.TileSpacing = 'compact';
    t2.Padding = 'compact';

    set(f6,'Position',[1 41 1536 840]) % full screen
    set(f7,'Position',[1 41 1536 840]) % full screen

    figs = [f1,f2,f3,f4,f5_arr,f6,f7];
end


function color_legend(groups,colors,loc)
    [~,first_idx] = ismember(categories(groups),groups);
    cols = colors(first_idx);
    h = gobjects([1 length(cols)]);
    for i=1:length(cols)
        h(i) = barh(NaN,NaN,cols{i});
    end
    set(h,{'DisplayName'},categories(groups));
    legend('location',loc)
end

function sensitivity_tornado_barh_inner(slope, param_names, colors, k)

    % set k=0 to show all parameters
    if k==0
        k = length(param_names);
    end

    % sort by highest sensitivity
    [~,sort_idx] = sort(abs(slope),'MissingPlacement','first');

    % only show the parameters with the k highest sensitivities.
    sort_idx_highest_k = sort_idx(end-k+1:end);
    slope_highest_k = slope(sort_idx_highest_k);
    colors_highest_k = colors(sort_idx_highest_k);
    params = categorical( param_names(sort_idx_highest_k), param_names(sort_idx_highest_k) );

    % bar chart
    for i=1:length(params)
        h = barh(params(i),slope_highest_k(i),colors_highest_k{i});
        h.Annotation.LegendInformation.IconDisplayStyle = 'off';
        hold on
    end
    xlim([-1.2 1.2] + fix(imrange(slope)))
    ax = gca;
    set(ax,'YGrid','on','XGrid','on')
    improvePlot

    % reduce fontsize
    all_axes = findobj(gcf().Children,'Type','Axes'); % all subplots in fig
    set(all_axes,'FontSize',12)
end

function f = sensitivity_tornado_barh(slope_LCOE, slope_Pvar, param_names, colors, groups, extra_title)
    
    % separate subplots for each objective
    f = figure;
    subplot 121
    sensitivity_tornado_barh_inner(slope_LCOE, param_names, colors, 25)
    title('At Minimum LCOE')
    
    subplot 122
    sensitivity_tornado_barh_inner(slope_Pvar, param_names, colors, 25)
    title('At Minimum Design Cost')
    sgtitle(['Normalized Sensitivities: ' extra_title])
    
    % legend
    color_legend(groups,colors,'southeast')
    set(f,"Position",[1.8 41.8 930 836.8])
    
    % both objectives on the same chart
%     figure
%     barh(categorical(param_names),[slope_LCOE, slope_Pvar])
%     legend('LCOE','P_{var}')
%     title(['Normalized Sensitivities' extra_title])
%     improvePlot
%     set(gca, 'FontSize', 12)
%     set(ax,'FontSize',12)
%     set(gca,'YGrid','on','XGrid','on')
end

function slope = get_slope(y, x)
    % y is a 2D (ie J*) or 3D (ie X*) array, and x is a 1D vector (ie ratios), 
    % then slope is similar to
    % (y(:,end,:) - y(:,1,:) ./ (x(end) - x(1))
    % but if y contains NaNs at those indices, then more inner indices are used.

    assert(size(x,2) == size(y,2))
    assert(ndims(y) <= 3)
    assert(size(x,1) == 1)
    x = x.'; % make row vector into col vector so indexing below works

    y_for_nans = y(:,:,1); % deal with case where ndims(y) > 2 - this works
    % because the presence of a NaN for 3D input (X*) is determined by the first two indices only

    % if there are no NaNs, idx_first will be all ones and idx_last will be
    % all size(result,2). Below calculates correct indices when there are NaNs.
    [~,sub_dim2_first] = max( ~isnan(y_for_nans), [], 2);
    flipped = flip(y_for_nans,2);
    [~,sub_dim2_last_flipped] = max( ~isnan(flipped), [], 2);
    sub_dim2_last = size(y,2) + 1 - sub_dim2_last_flipped;
    [sub_dim1, sub_dim3] = meshgrid(1:size(y,1), 1:size(y,3));
    sub_dim2_first = repmat(sub_dim2_first,[1 size(y,3)]);
    sub_dim2_last  = repmat(sub_dim2_last, [1 size(y,3)]);
    idx_first = sub2ind(size(y), sub_dim1.', sub_dim2_first, sub_dim3.');
    idx_last  = sub2ind(size(y), sub_dim1.', sub_dim2_last,  sub_dim3.');

    % using idx_first and idx_last, calculate slope
    y_first = y(idx_first);
    y_last = y(idx_last);
    x_first = x(sub_dim2_first);
    x_last = x(sub_dim2_last);
    slope = (y_last - y_first) ./ (x_last - x_first);
    
end

function [LCOE, X_LCOE, g_lambda_LCOE, ...
         P_var, X_Pvar, g_lambda_Pvar] = global_sens_reoptimize(x0, p, b, ...
                                                                    param_name, ratios, num_DVs, num_constr)
    % Brute force parameter sensitivity sweep (reoptimize for each param value)

    [LCOE,P_var] = deal(zeros(1,length(ratios)));
    [X_LCOE,X_Pvar] = deal(zeros(length(ratios),num_DVs));
    [g_lambda_LCOE,g_lambda_Pvar] = deal(zeros(length(ratios),num_constr));

    var_nom = p.(param_name);

    dry_run = false; % true means use random numbers, false means actual optimization

    parfor j=1:length(ratios)
        if ratios(j) ~=1   
            new_p = p;
            new_p.(param_name) = ratios(j) * var_nom;
            if dry_run
                [Xs_opt, obj_opt, flag] = deal(rand(num_DVs+1,2),rand(2,1),[1 1]);
                g_lambda_LCOE_tmp = rand(num_constr,1);
                g_lambda_Pvar_tmp = rand(num_constr,1);
            else
                [Xs_opt, obj_opt, flag, ~, lambdas] = gradient_optim(x0,new_p,b);
                g_lambda_LCOE_tmp = combine_g_lambda(lambdas(1),Xs_opt(:,1),new_p,b);
                g_lambda_Pvar_tmp = combine_g_lambda(lambdas(2),Xs_opt(:,2),new_p,b);
            end
            
            if flag(1) >= 1
                LCOE(j) = obj_opt(1);
                X_LCOE(j,:) = Xs_opt(1:end-1,1);
                g_lambda_LCOE(j,:) = g_lambda_LCOE_tmp;
            else
                [X_LCOE(j,:),LCOE(j),g_lambda_LCOE(j,:)] = deal(NaN);
            end
            if flag(2) >= 1
                P_var(j) = obj_opt(2);
                X_Pvar(j,:)= Xs_opt(1:end-1,2); 
                g_lambda_Pvar(j,:) = g_lambda_Pvar_tmp;
            else
                [X_Pvar(j,:), P_var(j),g_lambda_Pvar(j,:)] = deal(NaN);
            end
        end
    end   

end

function g_lambda = combine_g_lambda(lambda,X,p,b)
% returns a vector containing lambda where constraints are active,
% and g where constraints are inactive.

    all_lambda = [lambda.ineqnonlin; lambda.ineqlin; lambda.lower; lambda.upper];
    idx_g = all_lambda == 0;

    % sensitivities use g<=0 convention
    [~,~,g_nl] = simulation(X,p);
    [A_lin, b_lin] = lin_ineq_constraints(p);
    g_lin = A_lin*X(1:size(A_lin,2)) - b_lin;
    g_lb = b.X_mins - X(1:end-1);
    g_ub = X(1:end-1) - b.X_maxs;
    all_g = [-g_nl.'; g_lin; g_lb; g_ub];

    g_lambda = all_lambda;
    g_lambda(idx_g) = all_g(idx_g);

end

