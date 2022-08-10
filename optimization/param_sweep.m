% Brute force parameter sensitivity sweep (reoptimize for each param value)
%% Setup
clear;clc;%close all
dvar_names={'D_f','D_{s_{ratio}}', 'h_{f_{ratio}}','T_{s_{ratio}}', 'F_{max}','B_p','w_n','M'};

param_names = {'Hs','Hs_{struct}','T','T_{struct}','\sigma_y','\rho_m','E',...
            'cost_m','t_{ft}','t_{fr}','t_{fc}', 't_{fb}','t_{sr}','t_{dt}','D_{dt}','\theta_{dt}',...
            'FCR','N_{WEC}','eff_{pto}','D_d/D_s','T_s/D_s','h_d/D_s','T_f/h_f'};  % list of parameters to sweep
params = regexprep(param_names,'[{}\\]','');    % remove the curly braces and slashes
params = regexprep(params,'/','_over_');

groups = {'Dynamics','Structures',...
    'Economics','Geometry'};
param_groupings = [1 2 1 2 2 2 2 3 2 2 2 2 2 2 4 4 3 3 1 4 4 4 4];
color_groupings = {'b','y','g','r'};
colors = color_groupings(param_groupings);
cell2table([params',groups(param_groupings)'])

ratios = .8 : .1 : 1.2;
p = parameters();
b = var_bounds(p);
%%
% use the optimal x as x0 to speed up the sweeps
x0 = struct('D_f',b.D_f_nom,'D_s_ratio',b.D_s_ratio_nom,'h_f_ratio',...
        b.h_f_ratio_nom,'T_s_ratio',b.T_s_ratio_nom,'F_max',b.F_max_nom,...
        'D_int',b.B_p_nom,'w_n',b.w_n_nom,'M',b.M_nom);
x0_vec = gradient_optim(x0,p,b);
x0 = struct('D_f',x0_vec(1,1),'D_s_ratio',x0_vec(2,1),'h_f_ratio',x0_vec(3,1),...
    'T_s_ratio',x0_vec(4,1),'F_max',x0_vec(5,1),'D_int',x0_vec(6,1),'w_n',x0_vec(7,1));
x0(2) = struct('D_f',x0_vec(1,2),'D_s_ratio',x0_vec(2,2),'h_f_ratio',x0_vec(3,2),...
    'T_s_ratio',x0_vec(4,2),'F_max',x0_vec(5,2),'D_int',x0_vec(6,2),'w_n',x0_vec(7,2));

LCOE  = zeros(length(params),length(ratios));
P_var = zeros(length(params),length(ratios));
X = zeros(length(params), length(ratios), 8);
X_LCOE= zeros(length(params), length(ratios), 8);
X_Pvar= zeros(length(params), length(ratios), 8);

%% Run optimization
for i=1:length(params)
    p = parameters();
    var_nom = p.(params{i});
    for j=1:length(ratios)
        p.(params{i}) = ratios(j) * var_nom;
        [Xs_opt, obj_opt, flag] = gradient_optim(x0,p,b);
        % [Xs_opt, obj_opt, flag] = deal(rand(8,2),rand(2,1),[1 1]); %dry run
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
xlim([-2.5 1])
ax = gca;
set(ax,'YGrid','on')
title('LCOE')

subplot 122
for i=1:length(LCOE_params)
    barh(Pvar_params(i),slope_Pvar(i),colors{i})
    hold on 
end
set(gca,'YGrid','on')
title('c_v')
sgtitle('Normalized Sensitivities')
legend_idx = zeros(1,max(param_groupings));
for i=1:max(param_groupings)
    legend_idx(i) = find(param_groupings==i,1);
end
labels = repmat({''},size(param_groupings));
labels(legend_idx) = groups;
legend(labels)
improvePlot
xlim([-2.5 1])
set(gca, 'FontSize', 14)
set(ax,'FontSize',14)

% both objectives on the same chart
figure
barh(categorical(param_names),[slope_LCOE; slope_Pvar])
legend('LCOE','P_{var}')
title('Sensitivities')

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
