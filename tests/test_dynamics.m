classdef (SharedTestFixtures={ ...
        matlab.unittest.fixtures.CurrentFolderFixture('.')...
        }) test_dynamics < matlab.unittest.TestCase
    % class based unit tests, as in https://www.mathworks.com/help/matlab/matlab_prog/class-based-unit-tests.html
    
    properties (Constant)
        run_wecsim_tests = true;
    end

    properties
        errors_singlebody
        errors_multibody
        errors_report
        table
        figs
    end

    methods(TestClassSetup)
        % Shared setup for the entire test class

        function runNominalDynamics(testCase)
            if testCase.run_wecsim_tests
                warning('off','MATLAB:contour:ConstantData')
                t = tic;
                obj = Wecsim();
                obj = obj.run_all_from_load_if_possible();
                wecsim_runtime = toc(t);
                fprintf('WecSim took %g minutes',wecsim_runtime/60)
                warning('on','MATLAB:contour:ConstantData')

                err_structs = obj.end_result_struct.err_structs;
                testCase.errors_singlebody = err_structs{1};
                testCase.errors_multibody  = err_structs{2};
                testCase.errors_report     = err_structs{3};
                testCase.table             = obj.tab_array_display{1};
                testCase.figs              = obj.fig_array;
            end
        end
    end

    methods(TestMethodSetup)
        function check_whether_to_run(testCase)
            testCase.assumeTrue(testCase.run_wecsim_tests)
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
            table2latex(testCase.table,'test-results/table_13.tex')
        end

        function dynamicValidationFigures(testCase)
            for i = 1:length(testCase.figs)
                fig = testCase.figs(i);
                fig_name = ['Figure_WecSim_' num2str(i)];
                pdf_name = ['test-results/' fig_name];
                exportgraphics(fig,pdf_name)
                diagnostic = matlab.unittest.diagnostics.FigureDiagnostic(fig,'Prefix',[fig_name '_']);
                testCase.log(diagnostic);
            end
        end
    end

end