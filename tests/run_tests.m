import matlab.unittest.plugins.CodeCoveragePlugin
import matlab.unittest.plugins.codecoverage.CoverageReport

sourceCodeFolder = "mdocean";
addpath(genpath(sourceCodeFolder))
suite = testsuite("tests");
runner = testrunner("textoutput");

reportFolder = "coverageReport";
reportFormat = CoverageReport(reportFolder);

p = CodeCoveragePlugin.forFolder(sourceCodeFolder,"Producing",reportFormat,"IncludingSubfolders",true);
runner.addPlugin(p)
results = runner.run(suite);
open(fullfile("coverageReport","index.html"))