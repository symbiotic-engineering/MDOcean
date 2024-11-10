import matlab.unittest.plugins.CodeCoveragePlugin
import matlab.unittest.plugins.codecoverage.CoverageReport
import matlab.unittest.TestRunner;
import matlab.unittest.plugins.XMLPlugin;
import matlab.unittest.plugins.codecoverage.CoberturaFormat;
import matlab.unittest.plugins.TestReportPlugin

sourceCodeFolder = 'mdocean';
addpath(genpath(sourceCodeFolder))
suite = testsuite('tests');
runner = testrunner('textoutput');

date = datestr(now,'yyyy-mm-dd_HH.MM.SS');
cov_dir = ['code-coverage/' date];
test_dir = ['test-results/' date];
mkdir(cov_dir)
mkdir(test_dir)

reportFormat = [CoverageReport([cov_dir '/coverageReport']) CoberturaFormat([cov_dir '/coverage.xml'])];

p1 = CodeCoveragePlugin.forFolder({sourceCodeFolder}, 'IncludingSubfolders', true, 'Producing', reportFormat);
p2 = XMLPlugin.producingJUnitFormat([test_dir '/results.xml']);
p3 = TestReportPlugin.producingPDF([test_dir '/testreport.pdf']);

%runner.addPlugin(p1);
runner.addPlugin(p2);
runner.addPlugin(p3);

results = runner.runInParallel(suite);

%open([cov_dir '/coverageReport/index.html'])
open([test_dir '/testreport.pdf'])

display(results);
assertSuccess(results);
