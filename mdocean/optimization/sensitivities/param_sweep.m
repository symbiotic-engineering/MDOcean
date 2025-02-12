function [] = param_sweep(filename_uuid)

    %% Setup
    b = var_bounds();
    if nargin>0
        b.filename_uuid = filename_uuid;
    end
    [p,T] = parameters();
    
    param_names = T.name_pretty(T.sweep);  % list of parameters to sweep
    params = T.name(T.sweep);
    dvar_names = b.var_names_pretty(1:end-1);
    
    groups = categorical(T.subsystem(T.sweep));
    color_groupings = {'b','g','r','k','y'};
    colors = color_groupings(groups);
    
    
    %%
    % use the optimal x as x0 to speed up the sweeps
    % and obtain gradients
    x0 = b.X_start_struct;
    bypass = true; % read single obj optim results from file to speed up debugging
    if bypass
        load('gradient_optim_021225.mat')
    else
        [x0_vec, J0, ~, ~, lambdas, grads, hesses] = gradient_optim(x0,p,b);
    end
    
    num_DVs = length(dvar_names);
    num_constr = length(b.constraint_names);
    
    %% Obtain local sensitivity 
    disp('done optimizing, starting local param sensitivity')
    tic
    [par_x_star_par_p_local, ...
     dJstar_dp_local, dJdp_local, ...
     par_J_par_p_local, ...
     delta_p_change_activity_local] = local_sens_both_obj_all_param(x0_vec, p, params, ...
                                                                    lambdas, grads, hesses, num_constr);
    toc
    
    disp('done local, startng global param sensitivity')
    [par_x_star_par_p_global, dJstar_dp_global, ...
        delta_p_change_activity_global] = global_sens_all_param(params, param_names, x0_vec, dvar_names, J0, p, b, colors, groups);

    dJdp_combined = [par_J_par_p_local dJdp_local dJstar_dp_local dJstar_dp_global.'];
    dJdp_names = {'local partial','local total linear','local total quadratic','global'};

    % normalization
    p_val = T.value{T.sweep};
    dJdp_normalized = dJdp_combined.' .* p_val ./ J0(1);
    dxdp_local_normalized  = par_x_star_par_p_local.'  .* p_val ./ x0_vec(1:end-1,1);
    dxdp_global_normalized = par_x_star_par_p_global.' .* p_val ./ x0_vec(1:end-1,1);
    delta_p_local_normalized  = delta_p_change_activity_local  ./ p_val;
    delta_p_global_normalized = delta_p_change_activity_global ./ p_val;
    
    % color grid plots
    sensitivity_plot(dJdp_normalized, 'dJ/dp normalized', param_names, dJdp_names, ...
        'Parameters p', 'Type of Sensitivity')
    
    dxdp_plot(dxdp_local_normalized,  param_names, dvar_names)
    dxdp_plot(dxdp_global_normalized, param_names, dvar_names)
    
    delta_p_plot(delta_p_local_normalized,  b, dvar_names, param_names)
    delta_p_plot(delta_p_global_normalized, b, dvar_names, param_names)

    % bar chart plots
end

function delta_p_plot(delta_p_norm,b,dvar_names,param_names)
    constr_names  = [b.constraint_names_pretty b.lin_constraint_names_pretty dvar_names dvar_names];
    idx_slam = contains(constr_names,'Slamming');
    delta_p_norm_combine_slamming = delta_p_norm(:,~idx_slam);
    delta_p_norm_slam = delta_p_norm(:,idx_slam);
    [~,idx_delta_p_slam] = min(abs(delta_p_norm_slam),[],2);
    idx_delta_p_slam = sub2ind(size(delta_p_norm_slam),1:length(param_names),idx_delta_p_slam');
    delta_p_norm_combine_slamming(:,end+1) = delta_p_norm_slam(idx_delta_p_slam);
    
    titl = 'Normalized delta p to change constraint activity';
    xticks = [constr_names(~idx_slam),'Slamming'];
    yticks = param_names;
    xlab = 'Constraints';
    ylab = 'Parameters';
    sensitivity_plot(delta_p_norm_combine_slamming, titl, xticks, yticks, xlab, ylab)
    
end

function dxdp_plot(dxdp, param_names, dvar_names)
    titl = 'dx*/dp normalized';
    xticks = param_names;
    yticks = dvar_names;
    xlab = 'Parameters p';
    ylab = 'Design Variables x';
    sensitivity_plot(dxdp, titl, xticks, yticks, xlab, ylab)
end

function sensitivity_plot(matrix, titl, xticks, yticks, xlab, ylab)
    color_each_element(matrix)
    title(titl)
    set(gca,'XTickLabel',xticks)
    set(gca,'YTickLabel',yticks)
    clim([-10 10])
    colormap(bluewhitered)
    xlabel(xlab)
    ylabel(ylab)
end

%% Rerun optimization for global sensitivities
function [par_x_star_par_p_global, ...
          dJstar_dp_global, ...
          delta_p_change_activity_global] = global_sens_all_param(params, param_names, x0_vec, dvar_names, J0, p, b, colors, groups)
    %ratios = .8 : .1 : 1.2;
    ratios = [.9 .95 1 1.05 1.1];
    num_DVs = size(x0_vec,1)-1;
    
    [LCOE, P_var]    = deal(zeros(length(params), length(ratios)));
    [X_LCOE, X_Pvar] = deal(zeros(length(params), length(ratios), num_DVs));
    x0 = cell2struct(num2cell(x0_vec),b.var_names',1);

    for i=1:length(params)
        param_name = params{i};
    
        [LCOE(i,:),  X_LCOE(i,:,:), ...
         P_var(i,:), X_Pvar(i,:,:)] = global_sens_reoptimize(x0, p, b, param_name, ratios, num_DVs);
    
    end
    
    % fill in ratios == 1 (left blank since same as initial)
    X_LCOE_0_3d = permute(repmat(x0_vec(1:end-1,1).',length(params),1,sum(ratios==1)),[1 3 2]);
    X_Pvar_0_3d = permute(repmat(x0_vec(1:end-1,2).',length(params),1,sum(ratios==1)),[1 3 2]);
    LCOE(:, ratios==1) = J0(1);
    P_var(:,ratios==1) = J0(2);
    X_LCOE(:,ratios==1,:) = X_LCOE_0_3d;
    X_Pvar(:,ratios==1,:) = X_Pvar_0_3d;
    
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
    
    %% Plot each sensitivity
    col_nom = find(ratios==1);
    LCOE_nom = LCOE(1,col_nom);
    Pvar_nom = P_var(1,col_nom);
    X_LCOE_nom = X_LCOE(1,col_nom,:);
    X_Pvar_nom = X_Pvar(1,col_nom,:);
    X_LCOE_nom = X_LCOE_nom(:);
    X_Pvar_nom = X_Pvar_nom(:);
    
    figure(1)
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
    
    for i = 1:num_DVs
        figure(2)
        subplot(3,5,i)
        plot(ratios, X_LCOE(:,:,i)./X_LCOE_nom(i))
        xlabel ('Parameter ratio from nominal')
        ylabel ('X* ratio from nominal')
        title([dvar_names{i} ' - min LCOE'])
        %improvePlot
        grid on
        
        figure(3)
        subplot(3,5,i)
        plot(ratios, X_Pvar(:,:,i)./X_Pvar_nom(i))
        xlabel ('Parameter ratio from nominal')
        ylabel ('X* ratio from nominal')
        title([dvar_names{i} ' - min c_v'])
        %improvePlot
        grid on
    end
    legend(param_names)
    
    %% Tornado chart for overall slope - objective sensitivities
    
    slope_LCOE = get_slope(LCOE, ratios);
    slope_Pvar = get_slope(P_var, ratios);
    
    slope_LCOE_norm = slope_LCOE / LCOE_nom; % normalize
    slope_Pvar_norm = slope_Pvar / Pvar_nom;

    if all(~isfinite(slope_LCOE),'all') || all(~isfinite(slope_Pvar),'all')
        msg = ['All slopes are NaN, meaning all optimizations failed. ' ...
            'This might be because no feasible solution can be found. ' ...
            'slope_LCOE = ', num2str(slope_LCOE), ' and slope_Pvar = ', num2str(slope_Pvar)];
        error(msg)
    end
    
    sensitivity_tornado_barh(slope_LCOE_norm, slope_Pvar_norm, param_names, colors, groups)
    
    %% assign outputs - fixme ignoring slope_Pvar for now
    dJstar_dp_global = slope_LCOE;
    par_x_star_par_p_global = zeros(length(params),num_DVs); % fixme zero
    delta_p_change_activity_global = zeros(length(params),233+6+2*14);

    %% Tornado chart for overall slope - X* sensitivities
    
    % fixme: these slopes are incorrect, the dimensions are off!
    % slope_X_LCOE = get_slope(X_LCOE, ratios, X_LCOE_nom);
    % slope_X_Pvar = get_slope(X_Pvar, ratios, X_Pvar_nom);
    % 
    % %slope_X_LCOE = (X_LCOE(:,end) - X_LCOE(:,1))./X_LCOE_nom;
    % %slope_X_Pvar = (X_Pvar(:,end) - X_Pvar(:,1))./X_Pvar_nom;
    % 
    % % separate charts for each design variable, both objectives on same chart
    % figure
    % for i=1:num_DVs
    %     subplot(2,4,i)
    %     barh(categorical(param_names),[slope_X_LCOE(i,:);slope_X_Pvar(i,:)])
    %     title(dvar_names{i})
    %     improvePlot
    % end
    % legend('LCOE','c_v')
end


function sensitivity_tornado_barh(slope_LCOE, slope_Pvar, param_names, colors, groups)

    [~,LCOE_sort_idx] = sort(abs(slope_LCOE),'MissingPlacement','first');
    [~,Pvar_sort_idx] = sort(abs(slope_Pvar),'MissingPlacement','first');

    LCOE_params = categorical(param_names, param_names(LCOE_sort_idx));
    Pvar_params = categorical(param_names, param_names(Pvar_sort_idx));
    
    % separate charts for each objective
    figure
    subplot 121
    for i=1:length(LCOE_params)
        barh(LCOE_params(i),slope_LCOE(i),colors{i})
        hold on
    end
    xlim([-1 1] + fix(imrange(slope_LCOE)))
    ax = gca;
    set(ax,'YGrid','on','XGrid','on')
    title('LCOE')
    
    subplot 122
    for i=1:length(LCOE_params)
        barh(Pvar_params(i),slope_Pvar(i),colors{i})
        hold on 
    end
    xlim([-1 1] + fix(imrange(slope_Pvar)))
    set(gca,'YGrid','on','XGrid','on')
    title('c_v')
    sgtitle('Normalized Sensitivities')
    
    % legend
    [~,first_idx] = ismember(categories(groups),groups);
    labels = repmat({''},size(groups));         % create blank labels
    labels(first_idx) = categories(groups);     % leave most labels blank, except the first of each category
    legend(labels)
    
    improvePlot
    set(gca, 'FontSize', 14)
    set(ax,'FontSize',14)
    
    % both objectives on the same chart
    figure
    barh(categorical(param_names),[slope_LCOE; slope_Pvar])
    legend('LCOE','P_{var}')
    title('Sensitivities')
    improvePlot
    set(gca, 'FontSize', 12)
    set(ax,'FontSize',12)
    set(gca,'YGrid','on','XGrid','on')
end

function slope = get_slope(y_result, x)
    result_for_nans = y_result(:,:,1); % deal with case where ndims(result) > 2

    % if there are no NaNs, idx_first will be all ones and idx_last will be
    % all size(result,2). Below calculates correct indices when there are NaNs.
    [~,col_first] = max( ~isnan(result_for_nans), [], 2);
    flipped = flip(result_for_nans,2);
    [~,col_last_flipped] = max( ~isnan(flipped), [], 2);
    col_last = size(y_result,2) + 1 - col_last_flipped;

    idx_first = sub2ind(size(y_result), 1:size(y_result,1), col_first');
    idx_last  = sub2ind(size(y_result), 1:size(y_result,1), col_last');

    y_first = y_result(idx_first);
    y_last = y_result(idx_last);
    x_first = x(col_first);
    x_last = x(col_last);
    slope = (y_last - y_first) ./ (x_last - x_first);
    
end

function [LCOE, X_LCOE, P_var, X_Pvar] = global_sens_reoptimize(x0, p, b, param_name, ratios, num_DVs)
% Brute force parameter sensitivity sweep (reoptimize for each param value)

    [LCOE,P_var] = deal(zeros(1,length(ratios)));
    [X_LCOE,X_Pvar] = deal(zeros(length(ratios),num_DVs));

    var_nom = p.(param_name);

    for j=1:length(ratios)
        if ratios(j) ~=1   
            p.(param_name) = ratios(j) * var_nom;
            %[Xs_opt, obj_opt, flag] = gradient_optim(x0,p,b);
            [Xs_opt, obj_opt, flag] = deal(rand(num_DVs,2),rand(2,1),[1 1]); %dry run
            if flag(1) >= 1
                LCOE(j) = obj_opt(1);
                X_LCOE(j,:) = Xs_opt(:,1);
            else
                [X_LCOE(j,:),LCOE(j)] = deal(NaN);
            end
            if flag(2) >= 1
                P_var(j) = obj_opt(2);
                X_Pvar(j,:)= Xs_opt(:,2); 
            else
                [X_Pvar(j,:), P_var(j)] = deal(NaN);
            end
        end
    end   

end

