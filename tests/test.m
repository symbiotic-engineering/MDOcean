classdef test < matlab.unittest.TestCase
    properties
        feasible
        failed
        simulated
        actual
    end
    properties (TestParameter)
        field = fieldnames(validation_inputs());
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
    
    methods(Test)
        % Test methods
        
        function allFiguresRun(testCase)
            testCase.verifyReturnsTrue( @all_figures )
        end

        function validateNominalFeasible(testCase)
            testCase.onFailure( ['Nominal design violates these constraints: ', testCase.failed] );
            testCase.verifyTrue(testCase.feasible);
        end

        function validateNominal(testCase, field)
            sim = testCase.simulated.(field);
            act = testCase.actual.(field);
            testCase.verifyEqual(sim, act, 'RelTol',.1)
        end

        function validateNominalHydroCoeffs(testCase)
            % see dev/wamit_coeff_plot
        end

        function hydrodynamicLimitObeyed(testCase)
            % see lines 41-45 in dynamics.m
        end
    end
    
end