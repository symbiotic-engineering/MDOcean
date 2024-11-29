import matlab.unittest.plugins.CodeCoveragePlugin
import matlab.unittest.plugins.codecoverage.CoverageReport
import matlab.unittest.TestRunner;
import matlab.unittest.plugins.XMLPlugin;
import matlab.unittest.plugins.codecoverage.CoberturaFormat;
import matlab.unittest.plugins.TestReportPlugin
import matlab.unittest.selectors.HasName
import matlab.unittest.constraints.ContainsSubstring

sourceCodeFolder = 'mdocean';
addpath(genpath(sourceCodeFolder))

if exist('../WEC-Sim','dir')
    wecSimFolder = '../WEC-Sim/source';
    set_param(0, 'ErrorIfLoadNewModel', 'off')
    addpath(genpath(wecSimFolder))
end

suite = testsuite('tests');
suite = selectIf(suite,~HasName(ContainsSubstring("dynamics"))); % filter out wecsim tests

pwd
disp({suite.BaseFolder}.')
disp(path)
runner = testrunner('textoutput');

date = datestr(now,'yyyy-mm-dd_HH.MM.SS');
cov_dir = ['code-coverage/' date];
test_dir = ['test-results/' date];
mkdir(cov_dir)
mkdir(test_dir)

reportFormat = [CoverageReport([cov_dir '/coverageReport']) CoberturaFormat([cov_dir '/coverage.xml'])];
dirOut = dir(fullfile(sourceCodeFolder, '**', '*.m'));
codeFilePaths = string({dirOut.folder}) + filesep + string({dirOut.name});
% don't test coverage of generated files
filePathsToExclude = contains(codeFilePaths, 'generated');
codeFilePaths(filePathsToExclude) = [];

p1 = CodeCoveragePlugin.forFile(codeFilePaths, 'Producing', reportFormat);
p2 = XMLPlugin.producingJUnitFormat([test_dir '/junit.xml']);
p3 = TestReportPlugin.producingPDF([test_dir '/testreport.pdf']);

runner.addPlugin(p1);
runner.addPlugin(p2);
runner.addPlugin(p3);

results = runner.runInParallel(suite);

if ~batchStartupOptionUsed % don't open reports when running on CI server 
    open([cov_dir '/coverageReport/index.html'])
    open([test_dir '/testreport.pdf'])
end

display(results);
assertSuccess(results);
