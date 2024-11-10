classdef test < matlab.unittest.TestCase
    % class based unit tests, as in https://www.mathworks.com/help/matlab/matlab_prog/class-based-unit-tests.html
    properties
        feasible
        failed
        simulated
        actual
    end

    % inputs for tests, including passing tolerances
    properties (TestParameter)
        field = fieldnames(validation_inputs());
        rel_tol = {.1,.1,.1,.1,.25,.25,.25,.1,.1,.1,.1,.1};
        which_figs = test.enumerateFigs()
        which_tabs = test.enumerateTabs()
    end

    % helper methods to enumerate all figures and tables
    methods (Static)
        function which_fig_struct = enumerateFigs()
            [~,num_figs,num_tabs,fig_names,~] = all_figures( [],[] );
            which_figs = [1:num_figs zeros(1, num_tabs)];
            none = strcat(repmat({'none'},1,num_tabs), string(1:num_tabs));
            fig_names = matlab.lang.makeValidName([fig_names, none]);
            which_figs = num2cell(which_figs);
            which_fig_struct = cell2struct(which_figs,fig_names,2);
        end
        function which_tab_struct = enumerateTabs()
            [~,num_figs,num_tabs,~,tab_names] = all_figures( [],[] );
            which_tabs = [zeros(1,num_figs), 1:num_tabs];
            none = strcat(repmat({'none'},1,num_figs), string(1:num_figs));
            tab_names = matlab.lang.makeValidName([none, tab_names]);
            which_tabs = num2cell(which_tabs);
            which_tab_struct = cell2struct(which_tabs,tab_names,2);
        end
    end

    methods(TestClassSetup)
        % Shared setup for the entire test class

        function changeFolder(testCase)
            import matlab.unittest.fixtures.CurrentFolderFixture
            desiredFolder = '../mdocean';
            testCase.applyFixture(CurrentFolderFixture(desiredFolder))
        end

        function runNominalValidation(testCase)
            % this is a shared setup because the results are used by both
            % validateNominal and validateNominalFeasible
            [feas, fail, sim, act] = validate_nominal_RM3();
            testCase.feasible  = feas;
            testCase.failed    = fail;
            testCase.simulated = sim;
            testCase.actual    = act;
        end
    end
    
    % Test methods
    methods(Test, ParameterCombination='sequential')   
        function allFiguresRun(testCase, which_figs, which_tabs)
            success_criterion = all_figures(which_figs,which_tabs);
            if ~isempty(success_criterion)
                for i=1:length(success_criterion)
                    testCase.verifyGreaterThan(success_criterion{i},0);
                end
            end
        end

        function validateNominal(testCase, field, rel_tol)
            sim = testCase.simulated.(field);
            act = testCase.actual.(field);
            testCase.verifyEqual(sim, act, 'RelTol',rel_tol)
        end
    end

    methods(Test)
        function validateNominalFeasible(testCase)
            testCase.onFailure( ['Nominal design violates these constraints: ', testCase.failed] );
            testCase.verifyTrue(testCase.feasible);
        end

        function validateNominalHydroCoeffs(testCase)
            mean_err = hydro_coeff_err(false);
            testCase.verifyLessThanOrEqual(mean_err, 0.10)
        end

        function hydrodynamicLimitObeyed(testCase)
            ratio = check_max_CW();
            testCase.verifyLessThanOrEqual( ratio, 1 );
        end
    end
    
end