import matlab.unittest.plugins.CodeCoveragePlugin
import matlab.unittest.plugins.codecoverage.CoverageReport
import matlab.unittest.TestRunner;
import matlab.unittest.plugins.XMLPlugin;
import matlab.unittest.plugins.codecoverage.CoberturaFormat;

sourceCodeFolder = "mdocean";
addpath(genpath(sourceCodeFolder))
suite = testsuite("tests");
runner = testrunner("textoutput");

mkdir('code-coverage');
mkdir('test-results');

reportFolder = "code-coverage/coverageReport";
reportFormat = CoverageReport(reportFolder);

p = CodeCoveragePlugin.forFolder(sourceCodeFolder,"Producing",reportFormat,"IncludingSubfolders",true);
p = CodeCoveragePlugin.forFolder({'.'}, 'IncludingSubfolders', true, 'Producing', CoberturaFormat('code-coverage/coverage.xml'));
p2 = XMLPlugin.producingJUnitFormat('test-results/results.xml');

runner.addPlugin(p);
runner.addPlugin(p2);

results = runner.run(suite);
display(results);

open(fullfile("coverageReport","index.html"))

assertSuccess(results);
