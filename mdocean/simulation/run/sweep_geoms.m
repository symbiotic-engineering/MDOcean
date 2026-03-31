close all

n = 3; % number of points per dimension

p = parameters();

zero_to_one   = linspace(0.01,0.99,n);
zero_to_three = linspace(0.01,3,   n);

h     = [10 200]; % gives a range of m0h of 0.36-2.06 and 2.40-39.79 for default wave periods
a1_a2 = zero_to_one;
a2_h = zero_to_three;
d1_h  = zero_to_one;
d2_d1 = zero_to_one;
a3_a1 = 1 + zero_to_three;

[H, A1_A2, A2_H, ...
 D1_H, D2_D1, A3_A1] = ndgrid(h, a1_a2, a2_h, d1_h, d2_d1, a3_a1);

A1_H = A1_A2 .* A2_H;
D2_H = D2_D1 .* D1_H;
A3_H = A3_A1 .* A1_H;

A1 = A1_H .* H;
A2 = A2_H .* H;
A3 = A3_H .* H;
D1 = D1_H .* H;
D2 = D2_H .* H;

b = var_bounds();
X = [b.X_noms; 1];

rerun = false;

if rerun
    m0h_stored = zeros([length(p.T),size(A1)]);
    hydro_ratio_result = zeros([length(p.T),size(A1)]);

    for i=1:numel(A1)
        disp(['Running ' num2str(i) ' of ' num2str(numel(A1))])
    
        D_s = 2*A1(i);
        D_f = 2*A2(i);
        D_d = 2*A3(i);
        T_s = D1(i);
        T_f_2 = D2(i);
    
        % construct params and design vec
        p_i = p;
        p_i.h = H(i);
        p_i.T_s_over_D_s = T_s / D_s;
        p_i.D_d_over_D_s = D_d / D_s;
        X_i = X;
        X_i(strcmp(b.var_names,'D_s')) = D_s;
        X_i(strcmp(b.var_names,'D_f')) = D_f;
        X_i(strcmp(b.var_names,'h_s')) = T_s + 5;
        X_i(strcmp(b.var_names,'T_f_2')) = T_f_2;
        
        subs = {ind2sub(size(A1),i)};
        m0h_stored(:,subs{:}) = dispersion(2*pi./p.T, p_i.h, p.g) * p_i.h;
    
        % run simulation (drag and force/power sat are turned off inside check_max_CW)
        [hydro_ratio, ~, ...
         ~, ~, ~, ~, ~, out] = check_max_CW('', p_i, X_i, false);
        if i==1
            val = out;
        else
            val(i) = out;
        end
    
        hydro_ratio_result(:,subs{:}) = max(hydro_ratio);
    end
else
    load('sweep_geom_results.mat')
end


m0h_tmp = m0h_stored(1,:,:,:,:,:,:);
result_tmp = hydro_ratio_result(1,:,:,:,:,:,:);
result_disc = discretize(result_tmp,10,'categorical');
[result_disc_sorted, idx_sort] = sort(result_disc(:));
idx = idx_sort(~isundefined(result_disc_sorted));

figure
pp = parallelplot([H(idx),A1_A2(idx),A2_H(idx),D1_H(idx), D2_D1(idx), A3_A1(idx),m0h_tmp(idx),result_tmp(idx)],'GroupData',result_disc_sorted(~isundefined(result_disc_sorted)));
pp.CoordinateTickLabels = {'h','a_1/a_2','a_2/h','d_1/h','d_2/d_1','a_3/a_1','m_0h','CW/CW_{max}'};
pp.Color = parula(10);

%% scatter plot
figure
% use 0-1 vars for RBG
red = myresize(A1_A2, length(p.T));
green = myresize(D1_H, length(p.T));
blue = myresize(D2_D1, length(p.T));
color = [red, green, blue]; 
% use A3_A1 for size
size_var = myresize(A3_A1, length(p.T));
scatter(m0h_stored(:),hydro_ratio_result(:),size_var,color)
xlabel('m_0 h')
ylabel('CW/CW_{max}')
ylim([0 1])

%% line plot
num_m0h = length(p.T)*length(h);
[m0h_mat,order] = myreshape(m0h_stored,num_m0h);
result_mat = myreshape(hydro_ratio_result,num_m0h,order(:,1));

red2 = A1_A2(1,:,:,:,:,:);
green2 = D1_H(1,:,:,:,:,:);
blue2 = D2_D1(1,:,:,:,:,:);
color2 = [red2(:), green2(:), blue2(:)];
red_var_name = 'a_1/a_2';
green_var_name = 'd_1/h';
blue_var_name = 'd_2/d_1';

size_var2 = A3_A1(1,:,:,:,:,:);
size_mult = 3;
size_var_name = 'a_3/a_1';

marker_type_var = A2_H(1,:,:,:,:,:);
marker_var_name = 'a_2/h';
marker_types = {'o','x','v','s','+','.'};

fig = figure;
ax = gca();
num_lines = numel(red2);
h_data = gobjects(num_lines,1);
for i = 1:num_lines
    idx_marker = marker_type_var(i) == a2_h;
    marker = marker_types{idx_marker};
    h_data(i) = semilogx(ax, m0h_mat(:,i), result_mat(:,i), ['-' marker], 'MarkerSize', size_mult * size_var2(i));
    set(h_data(i).Annotation.LegendInformation,'IconDisplayStyle','off')
    hold on
end
ax.ColorOrder = color2;
xlabel('m_0 h')
ylabel('CW/CW_{max}')
ylim([0 1.5])
m0h_minmax = [min(m0h_mat(:)) max(m0h_mat(:))];
plot(ax,m0h_minmax,[1 1],'k--')

% numeric legends
h_marker = gobjects(length(a2_h),1);
for i=1:length(a2_h)
    h_marker(i) = plot(ax,NaN,NaN,['k' marker_types{i}],'DisplayName',num2str(a2_h(i)));
end
h_marker_leg = legend(h_marker,'AutoUpdate','off');
title(h_marker_leg,marker_var_name)

ah1 = axes('position',get(ax,'position'));
unique_size = unique(size_var2);
h_size = gobjects(length(unique_size),1);
for i=1:length(unique_size)
    h_size(i) = plot(ax,NaN,NaN,'ko','MarkerSize',size_mult*unique_size(i),'DisplayName',num2str(unique_size(i)));
    hold on
end
h_size_leg = legend(ah1,h_size);
title(h_size_leg,size_var_name)
ah1.Visible ='off';

% color legend
colorAxPos = [.67 .43];
mini_plot_size = [.2 .22];
colorAx = axes('Position',[colorAxPos mini_plot_size]);
box on
plotSVG(loadSVG('RGBCube_a.svg'));
set(colorAx,'Ydir','reverse')
set(colorAx,'XTickLabel',[],'YTickLabel',[],'XTick',[],'YTick',[])
red_pos = [-1.65, .75];
green_pos = [1.18, .75];
blue_pos = [.17, -1.2];
text(red_pos(1),red_pos(2),red_var_name)
text(green_pos(1),green_pos(2),green_var_name)
text(blue_pos(1),blue_pos(2),blue_var_name)
axis image
axis off
improvePlot

% main plot formatting that must happen after improvePlot
set([h_data; h_marker],'MarkerFaceColor','none')
fig.Position(3:4)  = [1000 600];
xlim(ax,m0h_minmax)

function [byfreq,order] = myreshape(A,len,order)
    reshaped = reshape(A,len,[]);

    if nargin<3
        [byfreq,order] = sort(reshaped);
    else
        byfreq = reshaped(order,:);
    end
end

function collapsed = myresize(A,len)
% Repeats A len times so size(shifted) = [len size(A)].
% Then collapses into a column vector.
% Used to add an extra dimension for the period onto geometric meshes.
    n_repeats = [ ones(1,ndims(A)), len ];
    repeated = repmat(A,n_repeats);
    shifted = shiftdim(repeated,ndims(A));
    collapsed = shifted(:);
end