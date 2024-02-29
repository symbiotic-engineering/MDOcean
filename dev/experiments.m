clear;clc;close all

p = parameters();
b = var_bounds();

% X = [ 20 10 30;     % outer diameter of float	
%       .3 .1 .5;     % inner diameter ratio of float	
%       1.5 1.2 2;     % outer diameter of reaction plate	
%       1 2 3;        % material	
%       10 1 20;      % Number of WECs in array	
%       10 5 50       % D_int	
%       2*pi/7 2*pi/8 2*pi/9
%         ];	

n = 20;
% X = [ b.D_f_nom, linspace(b.D_f_min, b.D_f_max, n);
%       b.D_s_ratio_nom, linspace(b.D_s_ratio_min, b.D_s_ratio_max, n);
%       b.D_or_ratio_nom, linspace(b.D_or_ratio_min, b.D_or_ratio_max, n);
%       1, ones([1 n]);
%       b.F_max_nom, linspace(b.F_max_min, b.F_max_max, n);
%       b.D_int_nom, linspace(b.D_int_min, b.D_int_max, n);
%       b.w_n_nom, linspace(b.w_n_min, b.w_n_max, n) ];
%   ratios = X./X(:,1);

ratios = logspace(log10(1/3),log10(3),n);
ratios = [1, ratios(ratios~=1)];
%X =  [b.D_f_nom, b.D_s_ratio_nom, b.D_or_ratio_nom, 1, b.F_max_nom, b.B_p_nom, b.w_n_nom]' * ratios;
X =  [b.D_f_nom b.D_s_ratio_nom b.h_f_ratio_nom b.T_s_ratio_nom b.F_max_nom b.B_p_nom b.w_n_nom]' * ratios;
X(4,:) = ones(1,n+1);

X_nom = X(:,1);	
design_size = size(X);	
var_num = design_size(1);	
var_depth = design_size(2);	
LCOE = X*inf;
P_var = X*inf;
opt_idx = zeros(var_num,1);	
recommended = zeros(var_num,2);	
number_runs = var_depth + (var_num-1)*(var_depth-1);	
failed = cell(number_runs,1);	
power = zeros(number_runs,1);	
FOS = zeros(number_runs,1);	
X_ins = zeros(number_runs, var_num);	
design = 0;	
for i = 1:var_num	
    X_in = X_nom;	
    for j = 1:var_depth	
        if i == 1 || j~=1	
            changed_entry = X(i,j);	
            if ~isnan(changed_entry)	
                design = design+1;	
                X_in(i) = changed_entry;	
                X_ins(design,:) = X_in;	
                [LCOE_temp, P_var_temp, B, FOS1Y, FOS2Y, FOS3Y, ...
            FOS_buckling, GM, power(design)] = simulation(X_in,p);	
                FOS(design) = min([FOS1Y,FOS2Y,FOS3Y,FOS_buckling]);
                [feasible, failed{design}] = is_feasible(B, FOS(design), GM, P_var_temp, p);	
                if feasible	
                    LCOE(i,j) = LCOE_temp;	
                    P_var(i,j) = P_var_temp;
                else	
                    LCOE(i,j) = NaN;
                    P_var(i,j) = NaN;
                end	
            end	
        end	
    end	
    [~, opt_idx(i)] = min(LCOE(i,:));	
    recommended(i,:) = [X(i,opt_idx(i)), opt_idx(i)];	
end	
[LCOE_op, ~, ~, B_op, FOS1Y_op,FOS2Y_op,FOS3Y_op,FOS_buckling_op,power_op] = simulation(recommended(:,1),p);	
[op_feasible, CC] = is_feasible(power_op, B_op, FOS1Y_op, power_op, p);	

% create table for display	
var_names = {'D_f',...    % outer diameter of float (m)	
            'D_s/D_f',... % inner diameter ratio of float (m)	
            'D_or/D_f',...      % outer diameter ratio of reaction plate (m)	
            'M',...         % material (-)	
            'F_max',...     % number of WECs in array (-)	
            'B_p',...     % internal damping of controller (Ns/m)	
            'w_n'};         % natural frequency (rad/s)
var_names_pretty = {'D_{f}',...    % outer diameter of float (m)	
            'D_s/D_{f}',... % inner diameter ratio of float (m)	
            'D_{or}/D_{f}',...      % outer diameter of reaction plate (m)	
            'M',...         % material (-)	
            'F_{max}',...     % number of WECs in array (-)	
            'D_{int}',...     % internal damping of controller (Ns/m)	
            '\omega_n'};         % natural frequency (rad/s)
results = array2table(X_ins, 'VariableNames', var_names);	
LCOE = LCOE';
P_var = P_var';
results = addvars(results, round(LCOE(LCOE~=Inf),1), round(P_var(P_var~=Inf)), failed, ...	
    'NewVariableNames', {'LCOE ($/kWh)','c_v (%)','ConstraintsFailed'});	
disp(results)

%% create pareto front
figure
plot(LCOE, P_var, '*')
xlabel('LCOE')
ylabel('P_{var}')
title('Design of Experiments Pareto Front')
legend(var_names_pretty)
improvePlot

%% sensitivities plot
[ratios_sorted,idx] = sort(ratios);
LCOE(1,:) = LCOE(1,1); % fill in nominal LCOE results for each DV where it wasn't repeatedly tested
P_var(1,:) = P_var(1,1);

figure
t = tiledlayout(2,1);
t.TileSpacing = 'compact';

ax1 = nexttile(1);
plot(ratios_sorted,LCOE(idx,:).')
ylabel('LCOE ($/kWh)')
grid on
ax2 = nexttile(2);
plot(ratios_sorted,P_var(idx,:).')
ylabel('Power c_v (%)')
grid on

title(t,'Design of Experiments Results','FontWeight','bold','FontSize',20)
l = legend(var_names_pretty);
l.Location = 'bestoutside';
%l.Position = [.75 .5 .15 .3];
grid on
linkaxes([ax1,ax2],'x');
xlabel('Design Variable Ratio (-)')
xticklabels(ax1,{})
xticks(ax2,xticks(ax1))
improvePlot