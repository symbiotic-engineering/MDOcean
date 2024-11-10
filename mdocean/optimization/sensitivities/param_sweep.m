function [] = param_sweep(filename_uuid)

% Brute force parameter sensitivity sweep (reoptimize for each param value)
%% Setup
clear;clc;%close all
b = var_bounds();
[p,T] = parameters();
b.filename_uuid = filename_uuid;

param_names = T.name_pretty(T.sweep);  % list of parameters to sweep
params = T.name(T.sweep);
dvar_names = b.var_names_pretty;

groups = categorical(T.subsystem(T.sweep));
color_groupings = {'b','g','r','k','y'};
colors = color_groupings(groups);

ratios = .8 : .1 : 1.2;
%%
% use the optimal x as x0 to speed up the sweeps
x0 = b.X_start_struct;
x0_vec = gradient_optim(x0,p,b);
x0 = cell2struct(num2cell(x0_vec),b.var_names',1);

num_DVs = length(dvar_names);
LCOE  = zeros(length(params), length(ratios));
P_var = zeros(length(params), length(ratios));
X_LCOE= zeros(length(params), length(ratios), num_DVs);
X_Pvar= zeros(length(params), length(ratios), num_DVs);

%% Run optimization
for i=1:length(params)
    p = parameters();
    p.filename_uuid = filename_uuid;
    var_nom = p.(params{i});
    for j=1:length(ratios)
        p.(params{i}) = ratios(j) * var_nom;
        [Xs_opt, obj_opt, flag] = gradient_optim(x0,p,b);
        %[Xs_opt, obj_opt, flag] = deal(rand(num_DVs,2),rand(2,1),[1 1]); %dry run
        if flag(1) >= 1
            LCOE(i,j) = obj_opt(1);
            X_LCOE(i,j,:) = Xs_opt(:,1);
        else
            [X_LCOE(i,j,:),LCOE(i,j)] = deal(NaN);
        end
        if flag(2) >= 1
            P_var(i,j) = obj_opt(2);
            X_Pvar(i,j,:)= Xs_opt(:,2); 
        else
            [X_Pvar(i,j,:), P_var(i,j)] = deal(NaN);
        end
    end   
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

for i = 1:7
    figure(2)
    subplot(2,4,i)
    plot(ratios, X_LCOE(:,:,i)./X_LCOE_nom(i))
    xlabel ('Parameter ratio from nominal')
    ylabel ('X* ratio from nominal')
    title([dvar_names{i} ' - min LCOE'])
    %improvePlot
    grid on
    
    figure(3)
    subplot(2,4,i)
    plot(ratios, X_Pvar(:,:,i)./X_Pvar_nom(i))
    xlabel ('Parameter ratio from nominal')
    ylabel ('X* ratio from nominal')
    title([dvar_names{i} ' - min c_v'])
    %improvePlot
    grid on
end
legend(param_names)

%% Tornado chart for overall slope - objective sensitivities

slope_LCOE = get_slope(LCOE, ratios, LCOE_nom);
slope_Pvar = get_slope(P_var, ratios, Pvar_nom);

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
% for i=1:7
%     subplot(2,4,i)
%     barh(categorical(param_names),[slope_X_LCOE(i,:);slope_X_Pvar(i,:)])
%     title(dvar_names{i})
%     improvePlot
% end
% legend('LCOE','c_v')
end

function slope = get_slope(y_result, x, y_nominal)
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
    slope = slope / y_nominal; % normalize
end
