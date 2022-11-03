classdef test < matlab.unittest.TestCase
    properties
        feasible
        failed
        simulated
        actual
    end

    properties (TestParameter)
        field = fieldnames(validation_inputs());
        rel_tol = {.1,.1,.1,.1,.25,.25,.25,.1,.1,.1,.1,.1};
        which_figs = test.enumerateFigs()
        which_tabs = test.enumerateTabs()
    end

    methods (Static)
        function which_figs = enumerateFigs()
            [num_figs, num_tabs] = all_figures( [],[] );
            which_figs = [1:num_figs zeros(1, num_tabs)];
            which_figs = num2cell(which_figs);
        end
        function which_tabs = enumerateTabs()
            [num_figs, num_tabs] = all_figures( [],[] );
            which_tabs = [zeros(1,num_figs), 1:num_tabs];
            which_tabs = num2cell(which_tabs);
        end
    end

    methods(TestClassSetup)
        % Shared setup for the entire test class

        function changeFolder(testCase)
            import matlab.unittest.fixtures.CurrentFolderFixture
            desiredFolder = '../mdocean';
            testCase.applyFixture(CurrentFolderFixture(desiredFolder))
        end

        function runValidation(testCase)
            [feas, fail, sim, act] = validate_nominal_RM3();
            testCase.feasible  = feas;
            testCase.failed    = fail;
            testCase.simulated = sim;
            testCase.actual    = act;
        end
    end
    
    % Test methods
    methods(Test, ParameterCombination='sequential')   
        function allFiguresRun(~, which_figs, which_tabs)
            all_figures(which_figs,which_tabs);
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