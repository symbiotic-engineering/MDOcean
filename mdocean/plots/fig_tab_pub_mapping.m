function [figs_in_RE, figs_in_AOR, tabs_in_RE, tabs_in_AOR] = fig_tab_pub_mapping()

%% RE figures

figs_in_RE = cell([1,num_figs_RE]);
figs_in_RE{1}  = 'ReadNonMatlabFigs.RM3_image';
% 2: modeling methodology
figs_in_RE{2}  = 'ReadNonMatlabFigs.methodology_overview';
figs_in_RE{3}  = 'ReadNonMatlabFigs.N2_diagram';
figs_in_RE{4}  = 'ReadNonMatlabFigs.dimensions';
% 3: optimization methodology
figs_in_RE{5} = 'ReadNonMatlabFigs.optimization_flowchart';
% 4: results
figs_in_RE{6} = 'DesignSpaceExploration.experiments';
figs_in_RE{7} = 'Comparison.overlaid_geometry';
figs_in_RE{8} = 'Comparison.overlaid_hydro_coeffs';
figs_in_RE{9} = 'Comparison.probability_CDF';
figs_in_RE{10} = 'GradientOptimFigFunc.single_obj_opt_power_matrix';
figs_in_RE{11} = 'GradientOptimFigFunc.lagrange_multipliers';
figs_in_RE{12} = 'GradientOptimFigFunc.dJ_dx_gradient';
figs_in_RE{13} = 'Multistart.multistart_convergence_tree';
figs_in_RE{14} = 'Multistart.multistart_bar chart';
figs_in_RE{15} = 'ParamSensitivities.re_optim_objective_tornado';
figs_in_RE{16} = 'ParamSensitivities.re_optim_design_tornado_J1';
figs_in_RE{17} = 'ParamSensitivities.re_optim_design_tornado_J2';
figs_in_RE{18} = 'ParetoFigFunc.pareto_front_with_design_images';
figs_in_RE{19} = 'ParetoFigFunc.pareto_front_LCOE_contours';
figs_in_RE{20} = 'ParetoFigFunc.pareto_heuristics';
figs_in_RE{21} = 'ParetoFigFunc.pareto_constraint_activity';
figs_in_RE{22} = 'AllFigCompare.runtime_bar_chart';
figs_in_RE{23} = 'ParetoSweep.sweep_num_seeds';
% appendix E - optimization process
figs_in_RE{24} = 'ParamSensitivities.post_optim_re_optim_objective_grid';
figs_in_RE{25} = 'ParamSensitivities.post_optim_design_grid';
figs_in_RE{26} = 'ParamSensitivities.re_optim_design_grid';
% appendix F - supplementary results
figs_in_RE{27} = 'GradientOptimFigFunc.normalized_gradient';
figs_in_RE{28} = 'GradientOptimFigFunc.single_obj_convergence';
figs_in_RE{29} = 'Multistart.multistart_parallel_coordinates';
% graphical abstract (unnumbered so at the end)
figs_in_RE{30}  = 'ReadNonMatlabFigs.graphical_abstract';

%% AOR figures

figs_in_AOR = cell([1,num_figs_AOR]);
figs_in_AOR{1}  = 'ReadNonMatlabFigs.RM3_image';
% 2: modeling methodology
figs_in_AOR{2}  = 'ReadNonMatlabFigs.methodology_overview';
figs_in_AOR{3}  = 'ReadNonMatlabFigs.N2_diagram';
figs_in_AOR{4}  = 'ReadNonMatlabFigs.dimensions';
figs_in_AOR{5}  = 'ReadNonMatlabFigs.MEEM_geometry';
figs_in_AOR{6}  = 'SparHydro.spar_added_mass';
figs_in_AOR{7}  = 'HydroCoeffFigFunc.hydro_coeff_err';
figs_in_AOR{8}  = 'DescFcns.saturation_desc_fcn';
figs_in_AOR{9}  = 'DescFcns.drag_desc_fcn';
figs_in_AOR{10} = 'RunSingleFigFunc.nominal_power_matrix';
figs_in_AOR{11} = 'Wecsim.WECSim_error_histograms_multibody';
figs_in_AOR{12} = 'ReadNonMatlabFigs.FBD';
figs_in_AOR{13} = 'Validation.cost_vs_N_WEC';
figs_in_AOR{14} = 'Runtime.sim_runtime';
figs_in_AOR{15} = 'Runtime.hydro_runtime';
figs_in_AOR{16} = 'Runtime.dynamics_runtime';
% appendix A- hydro
figs_in_AOR{17} = 'Meem.meem_regions';
figs_in_AOR{18} = 'Meem.meem_sparsity';
figs_in_AOR{19} = 'Meem.meem_validation';
figs_in_AOR{20} = 'Meem.meem_matching';
figs_in_AOR{21} = 'Meem.meem_convergence';
figs_in_AOR{22} = 'Meem.asymptotic_b_vector';
% appendix B - dynamics
figs_in_AOR{23} = 'ForceSaturationFigFunc.power_force_sensitivity';
figs_in_AOR{24} = 'RunSingleFigFunc.drag_convergence';
figs_in_AOR{25} = 'Slamming.slamming_amplitude';
figs_in_AOR{26} = 'RunSingleFigFunc.slamming_model_comparison';
figs_in_AOR{27} = 'Wecsim.wecsim_all_sea_states';
% appendix C - structures
figs_in_AOR{28} = 'ReadNonMatlabFigs.equivalent_stiffness';
figs_in_AOR{29} = 'ReadNonMatlabFigs.trapezoid';
figs_in_AOR{30} = 'ReadNonMatlabFigs.damping_plate_flowchart';
figs_in_AOR{31} = 'DampingPlateStructures.damping_plate_moment';
figs_in_AOR{32} = 'DampingPlateStructures.damping_plate_deflection';
figs_in_AOR{33} = 'DampingPlateStructures.damping_plate_aspect_ratio';
% appendix D - economics
% graphical abstract (unnumbered so at the end)
figs_in_AOR{34}  = 'ReadNonMatlabFigs.graphical_abstract';

%% TABLES

tabs_in_AOR = cell(1,num_tabs_AOR);
tabs_in_AOR{1} = 'Cost.cost_parameters';
tabs_in_AOR{2} = 'Validation.validation';

tabs_in_RE = cell(1,num_tabs_RE);
tabs_in_RE{1} = 'Constraints.constraints';
tabs_in_RE{2} = 'DesignVars.design_vars';
tabs_in_RE{3} = 'Parameters.parameters';
tabs_in_RE{4} = 'Comparison.optimal_design_vars';
tabs_in_RE{5} = 'Comparison.optimal_outputs';
tabs_in_RE{6} = 'LocationSensitivity.location_sensitivity';

end