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

mkdir('code-coverage');
mkdir('test-results');

reportFormat = [CoverageReport('code-coverage/coverageReport') CoberturaFormat('code-coverage/coverage.xml')];

p1 = CodeCoveragePlugin.forFolder({sourceCodeFolder}, 'IncludingSubfolders', true, 'Producing', reportFormat);
p2 = XMLPlugin.producingJUnitFormat('test-results/results.xml');
p3 = TestReportPlugin.producingPDF('test-results/testreport.pdf');

runner.addPlugin(p1);
runner.addPlugin(p2);
runner.addPlugin(p3);

results = runner.run(suite);

open('code-coverage/coverageReport/index.html')
open('test-results/testreport.pdf')

display(results);
assertSuccess(results);
