% Brute force parameter sensitivity sweep (reoptimize for each param value)
%% Setup
clear;clc;close all
dvar_names={'D_{f}','D_{s_{ratio}}', 'h_{f_{ratio}}','T_{s_{ratio}}', 'F_{max}','D_{int}','w_{n}','M'};
param_names = {'Hs','Hs_{struct}','T','T_{struct}','\sigma_y','\rho_m','E',...
            'cost_m','t_{fr}','t_{fc}', 't_{fb}','t_{sr}','B_{min}','FOS_{min}','D_{d_{min}}',...
            'FCR','N_{WEC}','D_d/D_s','T_s/D_s','h_d/D_s','T_f/h_f'};  % list of parameters to sweep
params = regexprep(param_names,'[{}\\]','');    % remove the curly braces and slashes
params = regexprep(params,'/','_over_');
%%
ratios = .8 : .1 : 1.2;
p = parameters();
b = var_bounds(p);

% use the optimal x as x0 to speed up the sweeps
x0 = struct('D_f',b.D_f_nom,'D_s_ratio',b.D_s_ratio_nom,'h_f_ratio',...
        b.h_f_ratio_nom,'T_s_ratio',b.T_s_ratio_nom,'F_max',b.F_max_nom,...
        'D_int',b.D_int_nom,'w_n',b.w_n_nom,'M',b.M_nom);
x0_vec = gradient_optim(x0,p,b);
x0 = struct('D_f',x0_vec(1),'D_s_ratio',x0_vec(2),'h_f_ratio',x0_vec(3),...
    'T_s_ratio',x0_vec(4),'F_max',x0_vec(5),'D_int',x0_vec(6),'w_n',x0_vec(7));
   
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
        if flag >= 1
            LCOE(i,j) = obj_opt(1);
            P_var(i,j) = obj_opt(2);
            X_LCOE(i,j,:) = Xs_opt(:,1);
            X_Pvar(i,j,:)= Xs_opt(:,2); 
        else
            [X_LCOE(i,j,:),X_Pvar(i,j,:),LCOE(i,j), P_var(i,j)] = deal(NaN);
        end
    end   
end
%% Plot each sensitivity
col_nom = find(ratios==1);
LCOE_nom = LCOE(1,col_nom);
Pvar_nom = P_var(1,col_nom);
X_LCOE_nom = X_LCOE(1,col_nom,:);
X_Pvar_nom = X_Pvar(1,col_nom,:);

figure
plot(ratios,LCOE/LCOE_nom)
xlabel('Parameter Ratio from nominal')
ylabel('LCOE ratio from nominal')
legend(param_names)
improvePlot
grid on

figure
plot(ratios,P_var/Pvar_nom)
xlabel('Parameter ratio from nominal')
ylabel('Power Variation ratio from nominal')
legend(param_names)
improvePlot
grid on
for i= 1:8
figure
plot(ratios, X_LCOE(:,:,i)./X_LCOE_nom(i))
xlabel ('Parameter ratio from nominal')
ylabel ('X* ratio from nominal')
title([dvar_names{i} ' - min LCOE'])
legend(param_names)
improvePlot
grid on

figure
plot(ratios, X_Pvar(:,:,i)./X_Pvar_nom(i))
xlabel ('Parameter ratio from nominal')
ylabel ('X* ratio from nominal')
title([dvar_names{i} ' - min c_v'])
legend(param_names)
improvePlot
grid on
end

%% Plot overall slope as tornado chart

slope_LCOE = (LCOE(:,end) - LCOE(:,1))./LCOE_nom;
slope_Pvar = (P_var(:,end) - P_var(:,1))./Pvar_nom;
slope_X_LCOE = (X_LCOE(:,end) - X_LCOE(:,1))./X_LCOE_nom;
slope_X_Pvar = (X_Pvar(:,end) - X_Pvar(:,1))./X_Pvar_nom;

% separate charts for each objective
figure
subplot 121
barh(categorical(param_names),slope_LCOE)
title('LCOE')
subplot 122
barh(categorical(param_names),slope_Pvar)
title('c_v')
sgtitle('Normalized Sensitivities')
improvePlot
 
%Assuming we want X* sensitivity plotted separately
for i=1:8
figure
barh(categorical(param_names),[slope_X_LCOE(:,i),slope_X_Pvar(:,i)])
title(dvar_names{i})
legend('LCOE','c_{v}')
improvePlot

% figure
% barh(categorical(var_names),slope_X_Pvar(:,i))
% title('X*, C_v Normalized Parameter Sensitivities')
% improvePlot


end
% both objectives on the same chart
figure
barh(categorical(param_names),[slope_LCOE slope_Pvar])
legend('LCOE','P_{var}')
title('Sensitivities')
