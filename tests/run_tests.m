import matlab.unittest.plugins.CodeCoveragePlugin
import matlab.unittest.plugins.codecoverage.CoverageReport
import matlab.unittest.TestRunner;
import matlab.unittest.plugins.XMLPlugin;
import matlab.unittest.plugins.codecoverage.CoberturaFormat;
import matlab.unittest.plugins.TestReportPlugin
import matlab.unittest.plugins.DiagnosticsRecordingPlugin
import matlab.unittest.selectors.HasName
import matlab.unittest.constraints.ContainsSubstring

sourceCodeFolder = 'mdocean';
addpath(genpath(sourceCodeFolder))

warning('off','MATLAB:nearlySingularMatrix')

run_wecsim_validation = false;

if run_wecsim_validation && exist('../WEC-Sim','dir')
    warning('off','MATLAB:contour:ConstantData')
    wecSimFolder = '../WEC-Sim/source';
    set_param(0, 'ErrorIfLoadNewModel', 'off')
    addpath(genpath(wecSimFolder))
    cd("mdocean")
    try
        tic
        [singlebody, multibody, report, tab] = validate_dynamics();
        wecsim_runtime = toc;
        fprintf('WecSim took %g minutes',wecsim_runtime/60)
        save('wecsimresults_all','singlebody','multibody','report','tab')
    catch err
        warning( "wecsim failed: " + newline + getReport(err,"extended"))
    end
    warning('on','MATLAB:contour:ConstantData')
    cd("..")
end

suite = testsuite('tests');
if ~run_wecsim_validation
    suite = selectIf(suite,~HasName(ContainsSubstring("dynamics"))); % filter out wecsim tests
end
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
p3 = TestReportPlugin.producingPDF([test_dir '/testreport.pdf'],'IncludingPassingDiagnostics',true);
p4 = DiagnosticsRecordingPlugin('IncludingPassingDiagnostics',true);

runner.addPlugin(p1);
runner.addPlugin(p2);
runner.addPlugin(p3);
runner.addPlugin(p4);

runner.ArtifactsRootFolder = test_dir;

tic
results = runner.runInParallel(suite);
clockTime = toc;

if ~batchStartupOptionUsed % don't open reports when running on CI server 
    open([cov_dir '/coverageReport/index.html'])
    open([test_dir '/testreport.pdf'])
end

fig = optim_time_bar_chart(suite,results);
save_pdf(fig,'test-results/Figure_30.pdf')

display(results);
fprintf('Testing took %g minutes',clockTime/60)
assertSuccess(results);
