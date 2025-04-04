classdef (SharedTestFixtures={ ...
        matlab.unittest.fixtures.CurrentFolderFixture('../mdocean')...
        }) test_dynamics < matlab.unittest.TestCase
    % class based unit tests, as in https://www.mathworks.com/help/matlab/matlab_prog/class-based-unit-tests.html
    
    properties
        errors_singlebody
        errors_multibody
        errors_report
        table
    end

    methods(TestClassSetup)
        % Shared setup for the entire test class

        function runNominalDynamics(testCase)
            %[singlebody, multibody, report, tab] = validate_dynamics();
            r = load('wecsimresults_all');
            testCase.errors_singlebody = r.singlebody;
            testCase.errors_multibody  = r.multibody;
            testCase.errors_report     = r.report;
            testCase.table             = r.tab;
        end
    end

    methods(Test)

        % Dynamic tests
        function validateSinglebodyWecsimBaseline(testCase)
            err = testCase.errors_singlebody.pct_error_baseline;
            testCase.verifyLessThanOrEqual(abs(err), 2);
            % match wecsim to within 2 percent under perfect assumptions
        end

        function validateSinglebodyWecsimTotal(testCase)
            err = testCase.errors_singlebody.pct_error_total;
            testCase.verifyLessThanOrEqual(abs(err), 10);
            % match wecsim to within 10 percent under assumptions used for optim
        end

        function validateMultibodyWecsimBaseline(testCase)
            err = testCase.errors_multibody.pct_error_baseline;
            testCase.verifyLessThanOrEqual(abs(err), 2);
            % match wecsim to within 2 percent under perfect assumptions
        end

        function validateMultibodyWecsimTotal(testCase)
            err = testCase.errors_multibody.pct_error_total;
            testCase.verifyLessThanOrEqual(abs(err), 10);
            % match wecsim to within 10 percent under assumptions used for optim
        end

        function validateMultibodyReportBaseline(testCase)
            err = testCase.errors_report.pct_error_baseline;
            testCase.verifyLessThanOrEqual(abs(err), 5);
            % match RM3 report to within 5 percent under perfect assumptions
        end

        function validateMultibodyReportTotal(testCase)
            err = testCase.errors_report.pct_error_total;
            testCase.verifyLessThanOrEqual(abs(err), 10);
            % match RM3 report to within 10 percent under assumptions used for optim
        end

        function dynamicValidationTable(testCase)
            diagnostic = matlab.unittest.diagnostics.DisplayDiagnostic(testCase.table);
            testCase.log(diagnostic);
        end
    end

end