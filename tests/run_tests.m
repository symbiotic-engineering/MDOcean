function run_tests(which_paper)
% run_tests  Run the MDOcean MATLAB test suite.
%
%   run_tests()           runs all tests (AOR + RE figures + dynamics)
%   run_tests('AOR')      runs AOR-paper figure/table tests plus all
%                         non-figure validation tests and test_dynamics
%   run_tests('RE')       runs only RE-paper figure/table tests

if nargin < 1
    which_paper = 'all';
end

import matlab.unittest.plugins.CodeCoveragePlugin
import matlab.unittest.plugins.codecoverage.CoverageReport
import matlab.unittest.TestRunner;
import matlab.unittest.plugins.XMLPlugin;
import matlab.unittest.plugins.codecoverage.CoberturaFormat;
import matlab.unittest.plugins.TestReportPlugin
import matlab.unittest.plugins.DiagnosticsRecordingPlugin
import matlab.unittest.plugins.LoggingPlugin
import matlab.unittest.selectors.HasName
import matlab.unittest.constraints.ContainsSubstring

add_mdocean_path();

if ~exist('drag_integral.mat', 'file')
    make_drag_integral
end
sourceCodeFolder = fullfile(fileparts(which('add_mdocean_path')), 'mdocean');

suite = testsuite('tests');

if ~strcmpi(which_paper, 'all')
    [figs, tabs] = get_fig_tab_names(which_paper, which_paper);
    names = [figs, tabs];
    % Build a selector matching any allFiguresRun test for this paper
    fig_selector = HasName(ContainsSubstring(names{1}));
    for i = 2:length(names)
        fig_selector = fig_selector | HasName(ContainsSubstring(names{i}));
    end
    if strcmpi(which_paper, 'AOR')
        % Include AOR figure tests + all non-allFiguresRun tests
        % (covers validateNominal* in test.m and all of test_dynamics)
        non_fig_selector = ~HasName(ContainsSubstring('allFiguresRun'));
        suite = suite.selectIf(fig_selector | non_fig_selector);
    else
        % RE: only the RE figure/table tests
        suite = suite.selectIf(fig_selector);
    end
end

runner = testrunner('textoutput');

if strcmpi(which_paper, 'all')
    base_test_dir = 'test-results';
    base_cov_dir = 'code-coverage';
else
    base_test_dir = ['test-results/' lower(which_paper)];
    base_cov_dir = ['code-coverage/' lower(which_paper)];
end
date = datestr(now,'yyyy-mm-dd_HH.MM.SS');
cov_dir = [base_cov_dir '/' date];
test_dir = [base_test_dir '/' date];
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
p3 = TestReportPlugin.producingPDF([test_dir '/testreport.pdf'],'IncludingPassingDiagnostics',true,'LoggingLevel',2);
p4 = DiagnosticsRecordingPlugin('IncludingPassingDiagnostics',true);
p5 = LoggingPlugin.withVerbosity(2);

runner.addPlugin(p1);
runner.addPlugin(p2);
runner.addPlugin(p3);
runner.addPlugin(p4);
runner.addPlugin(p5);

runner.ArtifactsRootFolder = test_dir;

t = tic;
results = runner.run(suite);    % parallel would make it slower because it 
                                % would run the test class setup on every worker
clockTime = toc(t);

if ~batchStartupOptionUsed % don't open reports when running on CI server 
    open([cov_dir '/coverageReport/index.html'])
    open([test_dir '/testreport.pdf'])
end

fig = unittest_time_bar_chart(suite,results);
exportgraphics(fig,[base_test_dir '/Figure_30_unittest.pdf'])

display(results);
fprintf('Testing took %g minutes',clockTime/60)
assertSuccess(results);
