clear;clc;close all

vars = {'F_max'}; % {'F_max','FCR','t_sft'};
ratios = [.1, 1:.5:8]; % .5 : .1 : 1.5
b = var_bounds(parameters());
x0 = struct('D_sft',b.D_sft_nom,'D_i_ratio',b.D_i_ratio_nom,'D_or',...
        b.D_or_nom,'N_WEC',b.N_WEC_nom,'D_int',b.D_int_nom,'w_n',b.w_n_nom);
   
LCOE  = zeros(length(vars),length(ratios));
P_var = zeros(length(vars),length(ratios));

for i=1:length(vars)
    p = parameters();
    var_nom = p.(vars{i});
    for j=1:length(ratios)
        p.(vars{i}) = ratios(j) * var_nom;
        [~, obj_opt, flag] = gradient_optim(x0,p,b);
        LCOE(i,j) = obj_opt(1);
        P_var(i,j) = obj_opt(2);
    end   
end
%%
figure
plot(ratios,LCOE)
xlabel('Ratio from nominal')
ylabel('LCOE')
legend(vars)
improvePlot
grid on

figure
plot(ratios,P_var)
xlabel('Ratio from nominal')
ylabel('Power Variation')
legend(vars)
improvePlot
grid on
%%
slope = LCOE(:,1) - LCOE(:,end)