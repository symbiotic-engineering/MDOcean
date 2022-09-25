classdef test < matlab.unittest.TestCase
    
    methods(TestClassSetup)
        % Shared setup for the entire test class
    end
    
    methods(TestMethodSetup)
        % Setup for each test
    end
    
    methods(Test)
        % Test methods
        
        function allFiguresRun(testCase)
            import matlab.unittest.fixtures.CurrentFolderFixture
            desiredFolder = '../mdocean';
            testCase.applyFixture(CurrentFolderFixture(desiredFolder))

            testCase.verifyReturnsTrue(@all_figures)
        end

        %function validation(testCase)
        %    testCase.verifyFail("Unimplemented test");
        %end
    end
    
end