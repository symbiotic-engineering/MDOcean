
clear
close all

[C_grid_data,percent_capacity_wecs_data,...
  cost_wec_thresh_data,x_grid_nonwec_data,... % outputs with WEC
  ...
  cost_wec_thresh_max_0,x_grid_nonwec_0,C_grid_0,... % outputs with no WEC
  ...
  LOC,ELEC,CARBON,...                       % grid inputs
  ...
  ZETA,OMEGA_RATIO,POWER_RATIO,COST_WEC] ... % design inputs
                                                = make_dummy_data();

sz = size(LOC);

delta_C_grid = C_grid_data - C_grid_0;
delta_cost_thresh = cost_wec_thresh_data - cost_wec_thresh_max_0;
delta_x_grid_nonwec = x_grid_nonwec_data - x_grid_nonwec_0;

% organize data
IN = [ZETA(:), OMEGA_RATIO(:), POWER_RATIO(:),COST_WEC(:)];
in_name = {'\zeta','\omega_n/\omega_p','P_{max}/P_{peak,unsat}','Cost WEC'};
OUT = [delta_C_grid(:), delta_cost_thresh(:), delta_x_grid_nonwec(:), percent_capacity_wecs_data(:)];
out_name = {'\DeltaCost Grid','\DeltaCost WEC Threshold','\DeltaCapacity Non-WEC','% Capacity WEC'};
GROUP = [LOC(:),ELEC(:),CARBON(:)];
group_name = {'Location','Electrification Scenario','Carbon Constraint'};
group_vals = {{'Northeast','California'},{'Ref Elec','Med Elec','High Elec'},{'Low CO_2','Med CO_2','High CO_2'}};

% group scatter plots
group_scatter_plot(IN, OUT, GROUP, in_name, out_name, group_name, group_vals)


% Attempt 1: beta indpendent of grid inputs (fit across all grid inputs)

% (step 1) beta 6,7,8: from fitting x_grid_nonwec_data - x_grid_nonwec_0
% as a function of ZETA, OMEGA_RATIO, and POWER_RATIO:
% x_grid_nonwec_model - x_grid_nonwec_0 = - beta(6)*ZETA - beta(7)*OMEGA_RATIO - beta(8)*POWER_RATIO;
in1 = [ZETA(:), OMEGA_RATIO(:), POWER_RATIO(:)];
out1 = delta_x_grid_nonwec(:);
ft1 = @(b,x) - b(1)*x(:,1) - b(2)*x(:,2) - b(3)*x(:,3);
guess1 = [1.5,1.5,1.5];
f1 = fitnlm(in1,out1,ft1,guess1,'CoefficientNames',{'b6','b7','b8'});

% (step 1 parallel) beta 3,4,5: from fitting cost_wec_thresh_data - cost_wec_thresh_max_0
% as a function of ZETA, OMEGA_RATIO, and POWER_RATIO:
% cost_wec_thresh_model - cost_wec_thresh_max_0 = - beta(3)*ZETA - beta(4)*OMEGA_RATIO - beta(5)*POWER_RATIO;
in2 = in1;
out2 = delta_cost_thresh(:);
ft2 = ft1;
guess2 = guess1;
f2 = fitnlm(in2,out2,ft2,guess2,'CoefficientNames',{'b3','b4','b5'});

% (step 2) beta 2: from fitting percent_capacity_wecs_data as a function of
% ratio_below_cost_thresh_fit, which is a function of beta3,4,5,
% cost_wec_thresh_max_0, ZETA, OMEGA_RATIO, and POWER_RATIO, COST_WEC:
% percent_capacity_wecs_model = min( beta(2) * ratio_below_cost_thresh_fit, 1);
cost_wec_thresh_fit = cost_wec_thresh_max_0 + reshape(f2.feval(in2),sz);
cost_wec_ratio_fit = cost_wec_thresh_fit ./ COST_WEC;
ratio_below_cost_thresh_fit = max( cost_wec_ratio_fit - 1, 0);

in3 = ratio_below_cost_thresh_fit(:);
out3 = percent_capacity_wecs_data(:);
ft3 = @(b,x) min( b * x, .99);
guess3 = 1.5;
f3 = fitnlm(in3,out3,ft3,guess3,'CoefficientNames',{'b2'});

% (step 3) beta 1: from fitting C_grid_data - C_grid_0 as a function of
% x_grid_wec_fit, which is a function of beta2,3,4,5,6,7,8,
% cost_wec_thresh_max_0, x_grid_nonwec_0, ZETA, OMEGA_RATIO, and POWER_RATIO, COST_WEC:
% C_grid_model - C_grid_0 = - beta(1) * x_grid_wec_fit;
percent_capacity_wecs_fit = reshape(f3.feval(in3),sz);
x_grid_nonwec_fit = x_grid_nonwec_0 + reshape(f1.feval(in1),sz);
x_grid_wec_fit = percent_capacity_wecs_fit .* x_grid_nonwec_fit ./ (1 - percent_capacity_wecs_fit);

in4 = x_grid_wec_fit(:);
out4 = delta_C_grid(:);
ft4 = @(b,x) -b * x;%'y~x1-1'
guess4 = guess3;
f4 = fitnlm(in4,out4,ft4,guess4,'CoefficientNames',{'b1'});

coeffs = [f4.Coefficients; f3.Coefficients; f2.Coefficients; f1.Coefficients]


function [C_grid_data,percent_capacity_wecs_data,...
          cost_wec_thresh_data,x_grid_nonwec_data,...
          cost_wec_thresh_max_0,x_grid_nonwec_0,C_grid_0,...
          LOC,ELEC,CARBON,...
          ZETA,OMEGA_RATIO,POWER_RATIO,COST_WEC] = make_dummy_data()

    location = [1,2];
    electrification_scenario = [1,2,3];
    carbon_constraint = [1,2,3];
    wave_cost = [1,2,3];
    zeta = [1,2,3];
    omega_n = [1,2,3];
    power_ratio = [.5 .8 1];
    [LOC,ELEC,CARBON,COST_WEC,ZETA,OMEGA_N,POWER_RATIO] = ndgrid(location,electrification_scenario,...
                                                    carbon_constraint,wave_cost,zeta,omega_n,power_ratio);

    num_beta = 8;
    num_fits = 4;
    beta = ones(1,num_beta); % fake underlying trend that fit should recover
    noise = rand([num_fits size(LOC)]);

    % constants dependent on non-wec-grid
    cost_wec_thresh_max_0 = .1*LOC + ELEC + CARBON;
    x_grid_nonwec_0 = LOC + ELEC + CARBON;
    C_grid_0 = LOC + ELEC + CARBON;
    OMEGA_P_0 = LOC;

    % model: what would come out of CEM if my model were perfect
    % data: what comes out of CEM (model + noise)
    % fit: what the fit to the data using the model predicts

    OMEGA_RATIO = OMEGA_N ./ OMEGA_P_0;
    cost_wec_thresh_model = cost_wec_thresh_max_0 - beta(3)*ZETA - beta(4)*OMEGA_RATIO - beta(5)*POWER_RATIO;
    cost_wec_ratio_model = cost_wec_thresh_model ./ COST_WEC;
    ratio_below_cost_thresh_model = max( cost_wec_ratio_model - 1, 0);
    percent_capacity_wecs_model = min( beta(2) * ratio_below_cost_thresh_model, .99);

    x_grid_nonwec_model = x_grid_nonwec_0 - beta(6)*ZETA - beta(7)*OMEGA_RATIO - beta(8)*POWER_RATIO;
    x_grid_wec_model = percent_capacity_wecs_model .* x_grid_nonwec_model ./ (1 - percent_capacity_wecs_model);
    C_grid_model = C_grid_0 - beta(1) * x_grid_wec_model;

    colons = repmat({':'},1,ndims(LOC));
    C_grid_data = C_grid_model + squeeze(noise(1,colons{:}));
    percent_capacity_wecs_data = percent_capacity_wecs_model + squeeze(noise(2,colons{:}));
    cost_wec_thresh_data = cost_wec_thresh_model + squeeze(noise(3,colons{:}));
    x_grid_nonwec_data = cost_wec_thresh_model + squeeze(noise(4,colons{:}));

end

function group_scatter_plot(IN,OUT,GROUP,in_name,out_name,group_name,group_vals)
    sym = 'osd*x+^v>'; % for grid scenario: electrification and carbon
    col = {'r','b'}; % for location
    col_cell = cellfun(@(x)repmat(x,[1 length(sym)]),col,'UniformOutput',false);
    color = strcat(col_cell{:});
    g = uigridlayout([1 2]);
    p  = uipanel(g);
    [h,ax,bigax] = gplotmatrix(p,IN,OUT,GROUP,color,sym,[],'off',[],in_name,out_name);
    p2 = uipanel(g);
    axLegend = axes(p2);
    for i=1:length(col)
        plot(axLegend,NaN,NaN,col{i},'DisplayName',group_vals{1}{i})
        hold(axLegend,'on')
    end
    for i=1:length(sym)
        num_carbons = length(group_vals{3});
        elec_name = group_vals{2}{ceil(i/num_carbons)};
        carbon_name = group_vals{3}{rem(i-1,num_carbons)+1};
        name = sprintf('Grid Scenario %i: %s, %s',i,elec_name,carbon_name);
        plot(axLegend,NaN,NaN,['k' sym(i)],'DisplayName',name)
    end
    legend(axLegend);
    axLegend.Visible = 'off';
end
