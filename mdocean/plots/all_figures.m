% Generate all figures used in the paper
function [fig_success,tab_success,...
          fig_output, tab_output,...
          fig_runtime,tab_runtime,...
          num_figs,num_tabs,fig_names,tab_names] = all_figures( which_figs, which_tabs, filename_uuid )

if nargin<3
    filename_uuid = ''; % required argument for anything running gradient_optim in parallel,
                        % since prob2struct needs unique filenames for code generation
end

date = datestr(now,'yyyy-mm-dd_HH.MM.SS');
save_folder = ['../test-results/' date '/'];
mkdir(save_folder)

num_figs = 36;
num_tabs = 8;
fig_names   = cell([1,num_figs]);
tab_names   = cell([1,num_tabs]);

if nargin==0
    % if run without arguments, show all figures and tables
    which_figs = 1:num_figs;
    which_tabs = 1:num_tabs;
end

fig_success = cell([1,length(which_figs)]);
tab_success = cell([1,length(which_tabs)]);

fig_output = gobjects(1, length(which_figs));
tab_output(1, 1:length(which_tabs)) = {table()};

fig_runtime = NaN([1,length(which_figs)]);
tab_runtime = NaN([1,length(which_tabs)]);

%% Define mapping from figures/tables to scripts
% figure_mapping = zeros(1,num_figs);
% figure_mapping(3:4) = 1; % sin desc function demo
% figure_mapping(5) = 2; % JPD
% figure_mapping(6:7) = 3; % pareto front
% figure_mapping(8) = 4; % parameter sensitivities
% figure_mapping(9:10) = 5; % overlaid comparison
% 
% table_mapping = zeros(1,num_tabs);
% table_mapping(1:2) = 6; % var bounds
% table_mapping(3) = 7; % parameters
% table_mapping(4) = 8; % validation
% table_mapping(5) = 5; % 3 design comparison
% table_mapping(6) = 9; % location sensitivity
% table_mapping(7) = 10; % convergence


%% non matlab figures
non_matlab_figs = [1:5 10:11 13:14];
non_matlab_fig_names = {'Fig. 1: RM3 image',...
                        'Fig. 2: Methodology overview',...
                        'Fig. 3: N2 diagram',...
                        'Fig. 4: Dimensions',...
                        'Fig. 5: MEEM geometry',...
                        'Fig. 10: WECSim error breakdown',...
                        'Fig. 11: FBD',...
                        'Fig. 13: sim runtime',...
                        'Fig. 14: Optimization flowchart'};
non_matlab_fig_files = {'geometry.png',...
                        'methods_flowchart_2_cropped.jpg',...
                        'xdsm.JPG',...
                        'dimensions.jpg',...
                        'MEEM-dims-basic-3.jpg',...
                        'Error_Accumulation_AEP.png',...
                        'structures_FBDs_4.jpg',...
                        'time_breakdown_3.png',...
                        'optim_process_flowchart_cropped.jpg'};

for i = 1:length(non_matlab_figs)
    fig_num = non_matlab_figs(i);
    fig_name = non_matlab_fig_names{i};
    fig_file = non_matlab_fig_files{i};

    fig_names{fig_num} = fig_name;
    if any(which_figs == fig_num)
        try
            t = tic;
            tmp_fig = figure;
            imshow(imread(fig_file),'Parent',axes(tmp_fig));
            tmp_fig.UserData = ['plots/non_matlab_figs/' fig_file(1:end-4) '.pdf'];
            fig_output(which_figs==fig_num) = tmp_fig;    
        catch err
            fig_success{which_figs==fig_num} = err;
        end
        fig_runtime(which_figs==fig_num) = toc(t);
    end
end

%% figure 6 - hydro coeffs vs freq
fig_names{6} = 'Fig. 6: hydro coeffs vs freq';
if any(which_figs == 6)
    try
        t = tic;
        hydro_coeff_err()
        fig6 = gcf;
        fig_output(which_figs==6) = fig6;
    catch err
        fig_success{which_figs==6} = err;
    end
    fig_runtime(which_figs==6) = toc(t);
end

%% figure 7, 8 - drag DF, saturation time signal
fig_names{7} = 'Fig. 7: drag describing function';
fig_names{8} = 'Fig. 8: force saturation time signal';
if any(which_figs == 7 | which_figs == 8)
    try
        t = tic;
        sin_desc_fcn_demo()
        fig7 = gcf;
        fig8 = figure(fig7.Number-2);
        fig_output(which_figs==7) = fig7;
        fig_output(which_figs==8) = fig8;
    catch err
        fig_success(which_figs == 7 | which_figs == 8) = {err};
    end
    fig_runtime(which_figs == 7 | which_figs == 8) = toc(t);
end

%% figure 9 - JPD multiplication
fig_names{9} = 'Fig. 9: JPD multiplication';
if any(which_figs == 9)
    try
        t = tic;
        p = parameters();
        b = var_bounds();
        X = [b.X_noms; 1];
        plot_power_matrix(X,p,b,filename_uuid)
        fig9 = gcf;
        fig_output(which_figs==9) = fig9;
    catch err
        fig_success{which_figs == 9} = err;
    end
    fig_runtime(which_figs==9) = toc(t);
end

%% figure 12 - cost vs N WEC
fig_names{12} = 'Fig. 12: cost vs N WEC';
if any(which_figs == 12) || any(which_tabs == 1)
    try
        t = tic;
        [~,~,~,~,tab1a,fig12] = validate_nominal_RM3('report');
        fig_output(which_figs==12) = fig12;
    catch err
        fig_success{which_figs == 12} = err;
    end
    fig_runtime(which_figs==12) = toc(t);
end

% %% figure 14-15 - wecsim validation histograms
% fig_names{14} = 'Fig. 14: WecSim histogram singlebody';
% fig_names{15} = 'Fig. 15: WecSim histogram multibody';
% if any(which_figs == 14 | which_figs == 15)
%     try
%         [~, ~, ~, ~, fig_singlebody, fig_multibody] = validate_dynamics();
%         fig_output(which_figs==14) = fig_singlebody;
%         fig_output(which_figs==15) = fig_multibody;
%     catch err
%         fig_success(which_figs == 14 | which_figs == 15) = {err};
%     end
% end

%% figure 15 - design space exploration
fig_names{15} = 'Fig. 15: design space exploration';
if any(which_figs == 15)
    try
        t = tic;
        experiments()
        fig15 = gcf;
        fig_output(which_figs==15) = fig15;
    catch err
        fig_success{which_figs == 15} = err;
    end
    fig_runtime(which_figs==15) = toc(t);
end

%% figure 16, 17, 18 - parameter sensitivities
fig_names{16} = 'Fig. 16: post-optimality and re-optimization dJ*/dp sensitivities';
fig_names{17} = 'Fig. 17: post-optimality dx*/dp sensitivities';
fig_names{18} = 'Fig. 18: re-optimization dx*/dp sensitivities';
if any(which_figs == 16 | which_figs == 17 | which_figs == 18)
    try
        t = tic;
        [runtimeLocal, runtimeGlobal] = param_sweep(filename_uuid);
        figTemp = gcf; % delta p global
        fig18 = figure(figTemp.Number - 2); % dx*/dp global
        fig17 = figure(figTemp.Number - 3); % dx*/dp local
        fig16 = figure(figTemp.Number - 5); % dJ*/dp combined
        fig_output(which_figs==16) = fig16;
        fig_output(which_figs==17) = fig17;
        fig_output(which_figs==18) = fig18;
        % fixme: 19 through 21 not implemented
        fig_runtime(which_figs==17) = runtimeLocal;
        fig_runtime(which_figs==18) = runtimeGlobal;
        fig_runtime(which_figs==16) = runtimeLocal + runtimeGlobal;
    catch err
        fig_success(which_figs == 16 | which_figs == 17 | which_figs == 18) = {err};
        fig_runtime(which_figs == 16 | which_figs == 17 | which_figs == 18) = toc(t);
    end
    
end

%% figure 22, 23, 24, 25 - pareto front, design heuristics
fig_names{22} = 'Fig. 22: pareto front with design images';
fig_names{23} = 'Fig. 23: design and objective heuristics';
fig_names{24} = 'Fig. 24: pareto front with LCOE contours';
fig_names{25} = 'Fig. 25: constraint activity';
if any(which_figs == 22 | which_figs == 23 | which_figs == 24 | which_figs == 25)
    try
        t = tic;
        pareto_search(filename_uuid);
        pareto_curve_heuristics()
        fig23 = gcf;
        fig22 = figure(fig23.Number - 2);
        fig24 = figure(fig23.Number - 3);
        fig25 = figure(fig23.Number - 5); % constraint activity
        fig25.Position = [1 41 1536 844.8000];
        fig_output(which_figs==22) = fig22;
        fig_output(which_figs==23) = fig23;
        fig_output(which_figs==24) = fig24;
        fig_output(which_figs==25) = fig25;
    catch err
        fig_success(which_figs == 22 | which_figs == 23 | which_figs == 24 | which_figs == 25) = {err};
    end
    fig_runtime(which_figs == 22 | which_figs == 23 | which_figs == 24 | which_figs == 25) = toc(t);
end

%% figure 26, 27, 28 - overlaid geometry, hydro coeffs, probability CDF
fig_names{26} = 'Fig. 26: overlaid geometry';
fig_names{27} = 'Fig. 27: overlaid hydro coeffs';
fig_names{28} = 'Fig. 28: probability CDF';
if any(which_figs == 26 | which_figs == 27 | which_figs == 28) || any(which_tabs == 5 | which_tabs == 6)
    try
        t = tic;
        [tab5,tab6] = compare(filename_uuid);
        n = gcf().Number;
        fig27 = figure(n);
        fig28 = figure(n-2);
        fig26 = figure(n-3);
        fig_output(which_figs==26) = fig26;
        fig_output(which_figs==27) = fig27;
        fig_output(which_figs==28) = fig28;
    catch err
        fig_success(which_figs == 26 | which_figs == 27 | which_figs == 28) = {err};
    end
    time = toc(t);
    fig_runtime(which_figs == 26 | which_figs == 27 | which_figs == 28) = time;
    tab_runtime(which_tabs == 5 | which_tabs == 6) = time;
end

%% figure 29 - fixme not implemented
fig_names{29} = 'Fig. 29: XXX';
if any(which_figs == 29)
    fig_output(which_figs==29) = figure;
end

%% figure 30 - asymptotic b vector
fig_names{30} = 'Fig. 30: asymptotic b vector';
if any(which_figs == 30)
    try
        t = tic;
        b_inf_numeric()
        fig30 = gcf;
        fig_output(which_figs==30) = fig30;
    catch err
        fig_success{which_figs==30} = err;
    end
    fig_runtime(which_figs==30) = toc(t);
end

%% figure 31-33 - fixme not implemented
if any( ismember(which_figs, 31:33) )
    fig_output(ismember(which_figs, 31:33)) = figure;
end

%% figure 34 to 36 - damping plate moment and deflection
fig_names{34} = 'Fig. 34: damping plate moment';
fig_names{35} = 'Fig. 35: damping plate deflection';
fig_names{36} = 'Fig. 36: damping plate plate aspect ratio';
if any(which_figs == 34 | which_figs == 35 | which_figs == 36)
    try
        t = tic;
        addpath('/damping-plate');
        BoedoPrantilAnnularPlate()
        fig36 = gcf;
        fig35 = figure(fig36.Number-1);
        fig34 = figure(fig36.Number-7);
        fig_output(which_figs==34) = fig34;
        fig_output(which_figs==35) = fig35;
        fig_output(which_figs==36) = fig36;
    catch err
        fig_success(which_figs == 34 | which_figs == 35 | which_figs == 36) = {err};
    end
    fig_runtime(which_figs == 34 | which_figs == 35 | which_figs == 36) = toc(t);
end

%% fixme: 37 to 45 not implemented, but still name the figures
empty_idx = cellfun(@isempty,fig_names);
empty_str = strcat('Fig._', string(find(empty_idx)));
[fig_names{empty_idx}] = deal(empty_str{:});
[fig_success{ismember(which_figs,find(empty_idx))}] = deal(MException('all_figures:not_implemented','Figure not implemented yet'));

%% table 12 - validation table
tab_names{1} = 'Tab. 12: validation against nominal';
if any(which_tabs == 1)
    try
        t = tic;
        % tab1a was generated with fig 12 above
        display(tab1a)
        [~,~,~,~,tab1b] = validate_nominal_RM3('wecsim');
        display(tab1b)
    
        % merge table 1a and 1b while preserving row order
        sharedCols = intersect(tab1a.Properties.VariableNames, tab1b.Properties.VariableNames);
        tab1a.RowNum = (1:length(tab1a.Properties.RowNames))';
        tab1b.RowNum = length(tab1a.Properties.RowNames) + (1:length(tab1b.Properties.RowNames))';
        
        tab1 = outerjoin(tab1a, tab1b, 'Keys', [{'RowNum'},sharedCols], 'MergeKeys', true);
        tab1 = removevars(tab1,'RowNum');
        tab1.Properties.RowNames = [tab1a.Properties.RowNames; tab1b.Properties.RowNames];
        display(tab1)
    
        tab_output{which_tabs==1} = tab1;
    
        vector_cols = {'capex','capex_design','J_capex_design','capex_struct','capex_PTO','opex','LCOE'};
        idx_remove = ismember(tab1.Properties.VariableNames,vector_cols); 
        tab1latex = rows2vars(tab1,'VariableNamingRule','preserve','DataVariables',~idx_remove);
        tab1latex.OriginalVariableNames = remove_underscores(modify_suffix(tab1latex.OriginalVariableNames));
        new_names = {'Variable','MDOcean','Actual','Error','MDOcean ','Actual ','Error '};
        tab1latex = renamevars(tab1latex, tab1latex.Properties.VariableNames, new_names);
        firstrow = '&\multicolumn{3}{c|}{DOE Report RM3 Design \cite{RM3}} & \multicolumn{3}{c}{WEC-Sim RM3 Design} \\';
        colspec = '>{\centering\arraybackslash}p{0.2\linewidth}|c|c|r|c|c|r';
        table2latex(tab1latex,[save_folder 'table_12.tex'],colspec,firstrow)
    catch err
        tab_success{which_tabs == 1} = err;
    end
    tab_runtime(which_tabs == 1) = toc(t);
end

%% table 15 - constraints table
tab_names{2} = 'Tab. 15: constraints';
if any(which_tabs == 2)
    try
        t = tic;
        b = var_bounds();
        tab2 = b.constraint_names';
        display(tab2)
        tab_output{which_tabs==2} = tab2;
    catch err
        tab_success{which_tabs == 2} = err;
    end
    tab_runtime(which_tabs == 2) = toc(t);
end

%% table 16 - design variables table
tab_names{3} = 'Tab. 16: design variables';
if any(which_tabs == 3)
    try
        t = tic;
        b = var_bounds();
        tab3 = array2table([b.X_mins b.X_noms b.X_maxs], ...
            'VariableNames',{'Mins','Noms','Maxs'}, 'RowNames', b.var_names(1:end-1));
        display(tab3)
        tab_output{which_tabs==3} = tab3;
    catch err
        tab_success{which_tabs == 3} = err;
    end
    tab_runtime(which_tabs == 3) = toc(t);
end

%% table 17 - parameters table
tab_names{4} = 'Tab. 17: parameters';
if any(which_tabs == 4)
    try
        t = tic;
        [~,tab4] = parameters();
        display(tab4)
        tab_output{which_tabs==4} = tab4;
    catch err
        tab_success{which_tabs == 4} = err;
    end
    tab_runtime(which_tabs == 4) = toc(t);
end

%% table 19 - optimal DVs for various designs
tab_names{5} = 'Tab. 19: optimal DVs for various designs';
if any(which_tabs == 5)
    try
        t = tic;
        % computation above with figures 27-29
        display(tab5);
        tab_output{which_tabs==5} = tab5;
        table2latex(tab5,[save_folder 'table_19.tex'])
    catch err
        tab_success{which_tabs == 5} = err;
    end
    tab_runtime(which_tabs == 5) = tab_runtime(which_tabs == 5) + toc(t);
end

%% table 20 - optimal outputs for various designs
tab_names{6} = 'Tab. 20: optimal outputs for various designs';
if any(which_tabs == 6)
    try
        t = tic;
        % computation above with figures 27-29
        display(tab6);
        tab_output{which_tabs==6} = tab6;
        table2latex(tab6,[save_folder 'table_20.tex'])
    catch err
        tab_success{which_tabs == 6} = err;
    end
    tab_runtime(which_tabs == 6) = tab_runtime(which_tabs == 6) + toc(t);
end

%% table 21 - convergence for different x0s
tab_names{7} = 'Tab. 21: convergence for different x0s';
if any(which_tabs == 7)
    try
        t = tic;
        tab7 = gradient_mult_x0(filename_uuid);
        tab_output{which_tabs==7} = tab7;
        table2latex(tab7,[save_folder 'table_21.tex'])
    catch err
        tab_success{which_tabs == 7} = err;
    end
    tab_runtime(which_tabs == 7) = toc(t);
end

%% table 22 - optimal DVs for 4 locations
tab_names{8} = 'Tab. 22: optimal DVs for 4 locations';
if any(which_tabs == 8)
    try
        t = tic;
        tab8 = location_sensitivity(filename_uuid);
        display(tab8);
        tab_output{which_tabs==8} = tab8;
        location_flags = str2double(tab8(strcmp(tab8.Row,'flag'),:).Variables);
        tab_success{which_tabs == 8} = {location_flags};
    
        idx_remove = ismember(tab8.Row,{'flag','Optimal Material index'});
        tab8latex = tab8(~idx_remove,:);
        colspec = ['>{\centering\arraybackslash}p{0.22\linewidth}' ...
                   '>{\centering\arraybackslash}p{0.08\linewidth}' ...
                   '>{\centering\arraybackslash}p{0.17\linewidth}' ...
                   '>{\centering\arraybackslash}p{0.17\linewidth}' ...
                   '>{\centering\arraybackslash}p{0.17\linewidth}' ...
                   '>{\centering\arraybackslash}p{0.18\linewidth}'];
        firstrow = '&& \multicolumn{4}{c}{Location}\\  \cline{3-6}';
        table2latex(tab8latex,[save_folder 'table_22.tex'],colspec,firstrow)
    catch err
        tab_success{which_tabs == 8} = err;
    end
    tab_runtime(which_tabs == 8) = toc(t);
end


end
