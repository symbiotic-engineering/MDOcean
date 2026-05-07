% Simple test for save_fig_with_diagnostic exportgraphics path
% Ensure MDOcean paths are on MATLAB path
repo_root = fileparts(fileparts(mfilename('fullpath')));
addpath(repo_root);
add_mdocean_path();

try
    mkdir('test-results');
catch
end

fig = figure('Visible','off');
plot(1:10, rand(1,10));
fig_name = 'test_exportgraphics';
pdf_prefix = 'test-results';

try
    diagnostic = save_fig_with_diagnostic(fig, fig_name, pdf_prefix);
    % expect a diagnostic object or DisplayDiagnostic
    if isa(diagnostic, 'matlab.unittest.diagnostics.Diagnostic') || isa(diagnostic, 'matlab.unittest.diagnostics.FigureDiagnostic') || isa(diagnostic, 'matlab.unittest.diagnostics.DisplayDiagnostic')
        fprintf('Diagnostic produced: %s\n', class(diagnostic));
    else
        fprintf('Returned type: %s\n', class(diagnostic));
    end
catch e
    disp(getReport(e));
    exit(1);
end

% Verify file exists
pdf_file = fullfile(pdf_prefix, [fig_name '.pdf']);
if exist(pdf_file, 'file') ~= 2
    fprintf('FAILED: PDF not created: %s\n', pdf_file);
    exit(2);
else
    fprintf('PASS: PDF created: %s\n', pdf_file);
    exit(0);
end
% Simple regression test for save_fig_with_diagnostic exportgraphics handling
% Reproduces the case where PDF doesn't exist and exportgraphics must be called
try
    fig = figure('Visible','off');
    fig_name = "damping_plate_flowchart";
    pdf_prefix = "test-results/";
    % ensure no accidental UserData value
    fig.UserData = [];
    diagnostic = save_fig_with_diagnostic(fig, fig_name, pdf_prefix);
    disp('save_fig_with_diagnostic: OK')
    close(fig);
catch ME
    if exist('fig','var') && isgraphics(fig)
        close(fig);
    end
    error('test_readnonmatlabfigs_exportgraphics:failed - %s', ME.message);
end
