clear;clc;close all

vars = {'F_max','FCR','t_sft'};
ratios = .5 : .1 : 1.5;
b = var_bounds(parameters());
x0 = struct('D_sft',b.D_sft_nom,'D_i_ratio',b.D_i_ratio_nom,'D_or_ratio',...
        b.D_or_ratio_nom,'N_WEC',b.N_WEC_nom,'D_int',b.D_int_nom,'w_n',b.w_n_nom);
   
LCOE = zeros(length(vars),length(ratios));
for i=1:length(vars)
    p = parameters();
    var_nom = p.(vars{i});
    for j=1:length(ratios)
        p.(vars{i}) = ratios(j) * var_nom;
        [X_opt, LCOE(i,j), flag] = gradient_optim(x0,p,b);
    end   
end
%%
figure
plot(ratios,LCOE)
xlabel('Ratio from nominal')
ylabel('LCOE')
legend({'F_{max}','FCR','t_{sft}'})
improvePlot
grid on
%%
slope = LCOE(:,1) - LCOE(:,end)