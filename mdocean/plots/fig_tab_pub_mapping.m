function [figs_in_RE, figs_in_AOR, tabs_in_RE, tabs_in_AOR] = fig_tab_pub_mapping()

    %% numbers
    num_figs_RE = 30;
    num_figs_AOR = 39;
    num_tabs_AOR = 2;
    num_tabs_RE = 6;

    %% RE figures

    figs_in_RE = cell([1, num_figs_RE]);
    figs_in_RE{1}  = 'ReadNonMatlabFigs.RM3_image';
    % 2: modeling methodology
    figs_in_RE{2}  = 'ReadNonMatlabFigs.methodology_overview';
    figs_in_RE{3}  = 'ReadNonMatlabFigs.N2_diagram';
    % archs
    figs_in_RE{4} = 'ReadNonMatlabFigs.optimization_flowchart';
    % 4: results
    figs_in_RE{5} = 'DesignSpaceExploration.experiments';
    figs_in_RE{6} = 'Comparison.overlaid_geometry';
    figs_in_RE{7} = 'Comparison.overlaid_hydro_coeffs';
    figs_in_RE{8} = 'Comparison.probability_CDF';
    figs_in_RE{9} = 'GradientOptimFigFunc.single_obj_opt_power_matrix';
    figs_in_RE{10} = 'GradientOptimFigFunc.lagrange_multipliers';
    figs_in_RE{11} = 'GradientOptimFigFunc.delta_x';
    figs_in_RE{12} = 'Multistart.multistart_convergence_tree';
    figs_in_RE{13} = 'Multistart.multistart_bar_chart';
    figs_in_RE{14} = 'ParamSensitivities.re_optim_objective_tornado';
    figs_in_RE{15} = 'ParamSensitivities.re_optim_design_tornado_J1';
    figs_in_RE{16} = 'ParamSensitivities.re_optim_design_tornado_J2';
    figs_in_RE{17} = 'ParetoFigFunc.pareto_front_with_design_images';
    figs_in_RE{18} = 'ParetoFigFunc.pareto_front_LCOE_contours';
    figs_in_RE{19} = 'ParetoFigFunc.pareto_heuristics';
    figs_in_RE{20} = 'ParetoFigFunc.pareto_constraint_activity';
    figs_in_RE{21} = 'ParetoFigFunc.pareto_damping_reactive';
    figs_in_RE{22} = 'AllFigCompare.runtime_bar_chart';
    figs_in_RE{23} = 'ParetoSweep.sweep_num_seeds';
    % appendix A - supplementary results
    figs_in_RE{24} = 'GradientOptimFigFunc.single_obj_convergence';
    figs_in_RE{25} = 'Multistart.multistart_parallel_coordinates';
    % appendix E - optimization process
    figs_in_RE{26} = 'ParamSensitivities.post_optim_re_optim_objective_grid';
    figs_in_RE{27} = 'ParamSensitivities.post_optim_design_grid';
    figs_in_RE{28} = 'ParamSensitivities.re_optim_design_grid';
    % graphical abstract (unnumbered so at the end)
    figs_in_RE{29}  = 'ReadNonMatlabFigs.graphical_abstract';

    %% AOR figures

    figs_in_AOR = cell([1, num_figs_AOR]);
    % 1: introduction
    figs_in_AOR{1}  = 'ReadNonMatlabFigs.RM3_image';
    figs_in_AOR{2}  = 'ReadNonMatlabFigs.methodology_overview';
    % 2: model structure
    figs_in_AOR{3}  = 'ReadNonMatlabFigs.N2_diagram';
    figs_in_AOR{4}  = 'ReadNonMatlabFigs.control_analysis_flowcharts';
    % 3: module details
    figs_in_AOR{5}  = 'ReadNonMatlabFigs.dimensions';
    figs_in_AOR{6}  = 'ReadNonMatlabFigs.MEEM_geometry';
    figs_in_AOR{7}  = 'SparHydro.spar_added_mass';
    figs_in_AOR{8}  = 'DescFcns.saturation_desc_fcn';
    figs_in_AOR{9}  = 'DescFcns.drag_desc_fcn';
    figs_in_AOR{10} = 'RunSingleFigFunc.nominal_power_matrix';
    figs_in_AOR{11} = 'ReadNonMatlabFigs.FBD';
    % 4: validation and benchmarking
    figs_in_AOR{12} = 'Wecsim.WECSim_error_histograms_multibody';
    figs_in_AOR{13} = 'HydroCoeffFigFunc.hydro_coeff_err';
    figs_in_AOR{14} = 'Validation.cost_vs_N_WEC';
    figs_in_AOR{15} = 'Runtime.sim_runtime';
    figs_in_AOR{16} = 'Runtime.hydro_runtime';
    figs_in_AOR{17} = 'Runtime.dynamics_runtime';
    % 5: insights and discussion
    figs_in_AOR{18} = 'DampingPlateStructures.damping_plate_aspect_ratio';
    figs_in_AOR{19} = 'ForceSaturationFigFunc.power_force_sensitivity';
    figs_in_AOR{20} = 'DesignSpaceExploration.experiments';
    figs_in_AOR{21} = 'RunSingleFigFunc.nominal_power_matrix'; % repeat of 10
    % appendix A- hydro
    figs_in_AOR{22} = 'Meem.meem_regions';
    figs_in_AOR{23} = 'Meem.meem_sparsity';
    figs_in_AOR{24} = 'Meem.meem_validation';
    figs_in_AOR{25} = 'Meem.meem_matching';
    figs_in_AOR{26} = 'Meem.meem_convergence';
    figs_in_AOR{27} = 'Meem.asymptotic_b_vector';
    % appendix B - dynamics
    figs_in_AOR{28} = 'RunSingleFigFunc.drag_convergence';
    figs_in_AOR{29} = 'Slamming.slamming_amplitude';
    figs_in_AOR{30} = 'RunSingleFigFunc.slamming_model_comparison';
    figs_in_AOR{31} = 'Wecsim.wecsim_all_sea_states';
    figs_in_AOR{31 + 6} = 'Wecsim.wecsim_all_sea_states_2';
    figs_in_AOR{31 + 7} = 'Wecsim.wecsim_all_sea_states_3';
    % appendix C - structures
    figs_in_AOR{32} = 'ReadNonMatlabFigs.equivalent_stiffness';
    figs_in_AOR{33} = 'ReadNonMatlabFigs.trapezoid';
    figs_in_AOR{34} = 'ReadNonMatlabFigs.damping_plate_flowchart';
    figs_in_AOR{35} = 'DampingPlateStructures.damping_plate_moment';
    figs_in_AOR{36} = 'DampingPlateStructures.damping_plate_deflection';
    % appendix D - economics
    % appendix E - parameters

    % graphical abstract (unnumbered so at the end)
    figs_in_AOR{39}  = 'ReadNonMatlabFigs.graphical_abstract';

    %% TABLES

    tabs_in_AOR = cell(1, num_tabs_AOR);
    tabs_in_AOR{1} = 'Cost.cost_parameters';
    tabs_in_AOR{2} = 'Validation.validation';

    tabs_in_RE = cell(1, num_tabs_RE);
    tabs_in_RE{1} = 'Constraints.constraints';
    tabs_in_RE{2} = 'DesignVars.design_vars';
    tabs_in_RE{3} = 'Parameters.parameters';
    tabs_in_RE{4} = 'Comparison.optimal_design_vars';
    tabs_in_RE{5} = 'Comparison.optimal_outputs';
    tabs_in_RE{6} = 'LocationSensitivity.location_sensitivity';

end
