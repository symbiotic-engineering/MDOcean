classdef (SharedTestFixtures={ ...
        matlab.unittest.fixtures.CurrentFolderFixture('../mdocean')}) ...
        test < matlab.unittest.TestCase
    % class based unit tests, as in https://www.mathworks.com/help/matlab/matlab_prog/class-based-unit-tests.html
    
    properties (Constant)
        run_slow_tests = true;
        slow_figs = [6 7 8];
        slow_tabs = 7;
    end

    properties
        feasible_report
        failed_report
        simulated_report
        actual_report

        feasible_wecsim
        failed_wecsim
        simulated_wecsim
        actual_wecsim

        uuid   
    end

    % inputs for tests, including passing tolerances
    properties (TestParameter)
        field_report = fieldnames(validation_inputs('report'));
        field_wecsim = fieldnames(validation_inputs('wecsim'));
        rel_tol_report = {.1,.1,.1,.1,.01,.01,.25,.25,.25,.1,.1,.1,.1,.1};
        rel_tol_wecsim = {.01,.01,.01,.01, 0.1,0.1};
        which_figs = test.enumerateFigs()
        which_tabs = test.enumerateTabs()
    end

    % helper methods to enumerate all figures and tables
    methods (Static)
        function which_fig_struct = enumerateFigs()
            [~,num_figs,num_tabs,fig_names,~] = all_figures( [],[] );

            if ~test.run_slow_tests
                num_tabs = num_tabs - length(test.slow_tabs);
            end

            which_figs_vec = [1:num_figs zeros(1, num_tabs)];
            none = strcat(repmat({'none'},1,num_tabs), string(1:num_tabs));
            fig_names = matlab.lang.makeValidName([fig_names, none]);

            if ~test.run_slow_tests
                idx_slow = ismember(which_figs_vec, test.slow_figs);
                which_figs_vec(idx_slow) = [];
                fig_names(idx_slow) = [];
            end

            which_figs_cell = num2cell(which_figs_vec);
            which_fig_struct = cell2struct(which_figs_cell,fig_names,2);
        end
        function which_tab_struct = enumerateTabs()
            [~,num_figs,num_tabs,~,tab_names] = all_figures( [],[] );

            if ~test.run_slow_tests
                num_figs = num_figs - length(test.slow_figs);
            end

            which_tabs_vec = [zeros(1,num_figs), 1:num_tabs];
            none = strcat(repmat({'none'},1,num_figs), string(1:num_figs));
            tab_names = matlab.lang.makeValidName([none, tab_names]);

            if ~test.run_slow_tests
                idx_slow = ismember(which_tabs_vec, test.slow_tabs);
                which_tabs_vec(idx_slow) = [];
                tab_names(idx_slow) = [];
            end

            which_tabs_cell = num2cell(which_tabs_vec);
            which_tab_struct = cell2struct(which_tabs_cell,tab_names,2);
        end
    end

    methods(TestClassSetup)
        % Shared setup for the entire test class

        function runNominalValidation(testCase)
            % this is a shared setup because the results are used by both
            % validateNominal and validateNominalFeasible
            [feas_r, fail_r, sim_r, act_r] = validate_nominal_RM3('report'); % report
            [feas_w, fail_w, sim_w, act_w] = validate_nominal_RM3('wecsim'); % wecsim
            
            testCase.feasible_report  = feas_r;
            testCase.failed_report    = fail_r;
            testCase.simulated_report = sim_r;
            testCase.actual_report    = act_r;

            testCase.feasible_wecsim  = feas_w;
            testCase.failed_wecsim    = fail_w;
            testCase.simulated_wecsim = sim_w;
            testCase.actual_wecsim    = act_w;
        end

        function generateUUID(testCase)
            % generate unique identifier for each parallel worker
            % required to prevent file overalps for generated code
            testCase.uuid = parallel.pool.Constant(@() char(matlab.lang.internal.uuid()));
        end
    end

    methods(TestClassTeardown)
        function deleteGeneratedFiles(testCase)
            % Create a wildcard pattern
            pattern = fullfile('**',['*' testCase.uuid.Value '*']);
            
            % Get a list of folders that match the pattern
            matchingFiles = dir(pattern);
            foldersToDelete = matchingFiles([matchingFiles.isdir]);
            foldersToDelete = fullfile({foldersToDelete.folder},{foldersToDelete.name});

            % Delete matching folders
            for d = 1:length(foldersToDelete)
                rmpath(foldersToDelete{d})
                rmdir(foldersToDelete{d},'s');
            end

        end
    end
    
    % Test methods
    methods(Test, ParameterCombination='sequential')   
        function allFiguresRun(testCase, which_figs, which_tabs)
            success_criterion = all_figures(which_figs,which_tabs,testCase.uuid.Value);
            if ~isempty(success_criterion)
                for i=1:length(success_criterion)
                    testCase.verifyGreaterThan(success_criterion{i},0);
                end
            end
        end

        function validateNominalReport(testCase, field_report, rel_tol_report)
            sim = testCase.simulated_report.(field_report);
            act = testCase.actual_report.(field_report);
            testCase.verifyEqual(sim, act, 'RelTol',rel_tol_report)
        end

        function validateNominalWecsim(testCase, field_wecsim, rel_tol_wecsim)
            sim = testCase.simulated_wecsim.(field_wecsim);
            act = testCase.actual_wecsim.(field_wecsim);
            testCase.verifyEqual(sim, act, 'RelTol',rel_tol_wecsim)
        end

    end

    methods(Test)
        % Static tests
        function validateNominalReportFeasible(testCase)
            testCase.onFailure( ['Nominal design violates these constraints: ', testCase.failed_report] );
            testCase.verifyTrue(testCase.feasible_report);
        end

        function validateNominalWecsimFeasible(testCase)
            testCase.onFailure( ['Nominal design violates these constraints: ', testCase.failed_wecsim] );
            testCase.verifyTrue(testCase.feasible_wecsim);
        end

        function validateNominalHydroCoeffs(testCase)
            mean_err = hydro_coeff_err(false);
            testCase.verifyLessThanOrEqual(mean_err, 0.10)
        end

        function hydrodynamicLimitObeyed(testCase)
            ratio = check_max_CW(testCase.uuid.Value);
            testCase.verifyLessThanOrEqual( ratio, 1 );
        end

    end
    
end