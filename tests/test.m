classdef test < matlab.unittest.TestCase
    properties
        feasible
        failed
        pct_error
    end

    methods(TestClassSetup)
        % Shared setup for the entire test class

        function changeFolder(testCase)
            import matlab.unittest.fixtures.CurrentFolderFixture
            desiredFolder = '../mdocean';
            testCase.applyFixture(CurrentFolderFixture(desiredFolder))
        end

        function runValidation(testCase)
            [feasible, failed, pct_error] = validate_nominal_RM3();
            testCase.feasible  = feasible;
            testCase.failed    = failed;
            testCase.pct_error = pct_error;
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

        function validateNominalMass(testCase)
            e = testCase.pct_error;
            mass_err = [e.mass_f e.mass_vc e.mass_rp e.mass_tot];
            testCase.verifyLessThanOrEqual(mass_err, .1)
        end

        function validateNominalCost(testCase)
            e = testCase.pct_error;
            cost_err = [e.capex e.opex];
            testCase.verifyLessThanOrEqual(cost_err, .05)
        end

        function validateNominalPower(testCase)
            e = testCase.pct_error;
            power_err = [e.power_avg e.power_max];
            testCase.verifyLessThanOrEqual(power_err, .1);
        end

        function validateNominalForce(testCase)
            e = testCase.pct_error;
            force_err = e.force_heave;
            testCase.verifyLessThanOrEqual(force_err, .1);
        end

        function validateNominalObjs(testCase)
            e = testCase.pct_error;
            obj_err = [e.LCOE e.c_v];
            testCase.verifyLessThanOrEqual(obj_err, .1);
        end

        function validateNominalHydroCoeffs(testCase)
            % see dev/wamit_coeff_plot
        end

        function hydrodynamicLimitObeyed(testCase)
            % see lines 41-45 in dynamics.m
        end
    end
    
end