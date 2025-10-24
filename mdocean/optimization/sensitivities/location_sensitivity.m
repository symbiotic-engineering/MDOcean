function [tab, h_pm, h_prob] = location_sensitivity(p,b)

files = {'Humboldt_California_Wave Resource _SAM CSV.csv',...
    'PacWave-North_Oregon_Wave-Resource.csv',...
    'PacWave-South_Oregon_Wave-Resource.csv',...
    'WETS_Hawaii_Wave-Resource.csv'};

depths = [45 50 65 55];
flux = [32.2 37 40.7 14.3];
storm = [NaN NaN NaN NaN];

X_opts = zeros(length(b.var_names),length(files));
obj_opts = zeros(1,length(files));
flags = zeros(1,length(files));
most_common_wave = cell(1,length(files));
BW = zeros(1,length(files));
BW_plot_on = false;

parfor i=1:length(files)
    new_p = p;
    
    jpd = trim_jpd(readmatrix(files{i}, 'Range', 'A3'));
    new_p.JPD = jpd(2:end,2:end);
    new_p.Hs = jpd(2:end,1);
    new_p.T = jpd(1,2:end);
    new_p.h = depths(i);

    new_b = fix_constraints(new_p,b);

    [~,idx_most_common] = max(new_p.JPD,[],'all');
    [row,col] = ind2sub(size(new_p.JPD),idx_most_common);
    most_common_wave(i) = {['$H_s = ' num2str(new_p.Hs(row)) '$m, $T_e=' num2str(new_p.T(col)) '$s']};
    BW(i) = round(find_BW(new_p.Hs,new_p.T,new_p.JPD,BW_plot_on),2);
  
    X = new_b.X_start_struct;
    
    which_obj = 1; % only optimize LCOE
    [X_opts(:,i), obj_opts(i), flags(i)]  = gradient_optim(X,new_p,new_b,which_obj);
    
    h_pm = plot_power_matrix(X_opts(:,i),new_p,b,b.filename_uuid);
    h_prob = figure(h_pm.Number + 1);
    power_PDF(X_opts(:,i),new_p)
    hold on
end
figure(h_pm.Number + 1)
locs = {'Humboldt, CA','PacWave North, OR', 'PacWave South, OR','Wave Energy Test Site, HI'};
legend(locs)

site_info = [flux;BW;storm;depths];
site_info_vars = {'$J$','$BW$','$H_{s,storm},T_{e,storm}$','$h$'};
site_info_descs = {'Incident energy flux (kW/m)',...
    'Half-power bandwidth (rad/s)','Storm sea states (m, s)','Water depth (m)'};

design_var_names = cellfun(@(x) ['$' x, '^*$'], b.var_names_pretty, 'UniformOutput', false);
design_var_descs = cellfun(@(x) ['Opt. ' x], lower(b.var_descs), 'UniformOutput', false);

symbols = [site_info_vars,  design_var_names,{'$LCOE^*$','flag'}];
row_descs = [site_info_descs, design_var_descs,{'Optimal levelized cost of energy (\\$/kWh)','flag'}];

first_row = array2table(string([{'-'},most_common_wave]),'VariableNames',[{'Symbol'},locs],'RowNames',{'Most frequent wave'});

start_tab = array2table([site_info;X_opts;obj_opts;flags],'VariableNames',locs,'RowNames',row_descs);
symbol_col = array2table(string(symbols).','VariableNames',{'Symbol'},'RowNames',row_descs);

main_tab =  horzcat(symbol_col,start_tab);
tab = vertcat(first_row,main_tab);

%% try Hawaii with California design
X_cali = X_opts(:,1);

LCOE_hawaii_with_cali_design = simulation(X_cali,p)
pct_diff = (LCOE_hawaii_with_cali_design - obj_opts(4)) / obj_opts(4)

end

function delta_w = find_BW(Hs,Te,JPD,plotOn)
    [T_mesh,Hs_mesh] = meshgrid(Te,Hs);
    energy_weighted_JPD = JPD .* T_mesh .* Hs_mesh.^2;
    energy_weighted_sum = sum(energy_weighted_JPD,1);
    w = 2*pi./Te;

    % flip both so w increases with idx
    w = fliplr(w);
    P = fliplr(energy_weighted_sum);

    % full width at half max
    [max_power,idx_max] = max(P);
    P_half = max_power/2;
    w_below_pk = w(1:idx_max);
    P_below_pk = P(1:idx_max);
    w_above_pk = w(idx_max:end);
    P_above_pk = P(idx_max:end);
    w_half_below = interp1(P_below_pk(P_below_pk>0),w_below_pk(P_below_pk>0),P_half);
    w_half_above = interp1(P_above_pk(P_above_pk>0),w_above_pk(P_above_pk>0),P_half);
    delta_w = w_half_above - w_half_below;

    if plotOn
        figure
        plot(w,P)
        hold on
        plot(w(idx_max),max_power,'mp')
        plot(w_half_below,P_half,'ro')
        plot(w_half_above,P_half,'go')
    end

end

function b = fix_constraints(p,b)
    num_non_sea_state_constraints = sum(~contains(b.constraint_names,'slamming'));
    desired_length = num_non_sea_state_constraints + numel(p.JPD) + length(p.T_struct);
    len = length(b.constraint_names);
    if len < desired_length
        b.constraint_names(end+1 : end+desired_length-len) = {'slamming'};
        b.constraint_names_pretty = remove_underscores(b.constraint_names);
    elseif len > desired_length
        b.constraint_names = b.constraint_names(1:desired_length);
        b.constraint_names_pretty = remove_underscores(b.constraint_names);
    end
end
