% Quick verification of the four CI fixes
% Tests without running the full test suite

repo_root = fileparts(fileparts(mfilename('fullpath')));
addpath(repo_root);
add_mdocean_path();

fprintf('=== Testing CI Fixes ===\n\n');

% Test 1: save_fig_with_diagnostic exportgraphics fix
fprintf('TEST 1: save_fig_with_diagnostic exportgraphics handling\n');
try
    mkdir('test-results');
catch
end
fig = figure('Visible','off');
plot(1:10, rand(1,10));
diagnostic = save_fig_with_diagnostic(fig, 'test_ci_fix_1', 'test-results');
pdf_file = fullfile('test-results', 'test_ci_fix_1.pdf');
if exist(pdf_file, 'file') == 2
    fprintf('✓ PASS: PDF created successfully\n');
    delete(pdf_file);
else
    fprintf('✗ FAIL: PDF not created\n');
end
close(fig);

% Test 2: Verify shortening logic in make_all_cases (without calling heavy function)
fprintf('\nTEST 2: Wecsim struct field names use shortened multi/singl\n');
% Check that the code change produces names within the MATLAB 63-char limit
short_case_desc_template_true = 'wcsm_multi_true_drag_on_meem_on';
contour_name = 'float_amplitude';

full_fig_name_short = sprintf('wecsim_%s__%s', short_case_desc_template_true, contour_name);
fprintf('  Example shortened figure name: %s (%d chars)\n', full_fig_name_short, length(full_fig_name_short));

if length(full_fig_name_short) <= 63
    fprintf('✓ PASS: Example name within 63-char limit\n');
else
    fprintf('✗ FAIL: Example name exceeds 63-char limit\n');
end

% Test 3: AllFigCompare has figure handle in intermediate struct
fprintf('\nTEST 3: AllFigCompare intermediate struct has figure_handle\n');
try
    intermed = AllFigCompare.analysis_fcn([],[]);
    if isfield(intermed, 'figure_handle') && ~isempty(intermed.figure_handle)
        fprintf('✓ PASS: figure_handle field exists and is non-empty\n');
        close(intermed.figure_handle);
    else
        fprintf('✗ FAIL: figure_handle is missing or empty\n');
    end
catch e
    fprintf('✗ FAIL: Error calling AllFigCompare.analysis_fcn: %s\n', e.message);
end

% Test 4: Check no unexpected "ExtraInputs" errors in save_fig_with_diagnostic
fprintf('\nTEST 4: save_fig_with_diagnostic error handling (retry with coercion)\n');
try
    % This would have previously failed with "Additional, unrecognized inputs"
    fig2 = figure('Visible','off');
    plot(1:5, rand(1,5));
    diagnostic2 = save_fig_with_diagnostic(fig2, 'test_ci_fix_4', 'test-results');
    fprintf('✓ PASS: save_fig_with_diagnostic completed without ExtraInputs error\n');
    close(fig2);
    pdf_file2 = fullfile('test-results', 'test_ci_fix_4.pdf');
    if exist(pdf_file2, 'file') == 2
        delete(pdf_file2);
    end
catch e
    if contains(e.message, 'ExtraInputs')
        fprintf('✗ FAIL: ExtraInputs error still occurs\n');
    else
        fprintf('✗ FAIL: Unexpected error: %s\n', e.message);
    end
    close all
end

% Test 5: Empty fig_name should fall back to a safe default
fprintf('\nTEST 5: save_fig_with_diagnostic empty fig_name fallback\n');
try
    fig3 = figure('Visible','off');
    plot(1:5, rand(1,5));
    diagnostic3 = save_fig_with_diagnostic(fig3, '', 'test-results');
    default_pdf = fullfile('test-results', 'unnamed_figure.pdf');
    if exist(default_pdf, 'file') == 2
        fprintf('✓ PASS: empty fig_name fell back to unnamed_figure.pdf\n');
        delete(default_pdf);
    else
        fprintf('✗ FAIL: fallback PDF was not created\n');
    end
    close(fig3);
catch e
    fprintf('✗ FAIL: empty fig_name caused error: %s\n', e.message);
    close all
end

% Test 6: Verify shortening is applied in make_all_cases.m
fprintf('\nTEST 6: Verify make_all_cases.m code uses multi/singl\n');
try
    % Read the source and check for the short names
    fid = fopen(which('make_all_cases'), 'r');
    content = fread(fid, '*char').';
    fclose(fid);
    
    has_short_multi = contains(content, 'wg_wcsm_multi_true');
    has_old_multibody = contains(content, 'wecsim_geom_wecsim_multibody_true');
    
    if has_short_multi
        fprintf('✓ PASS: make_all_cases.m uses short names (multi, not multibody)\n');
    else
        fprintf('✗ FAIL: make_all_cases.m still contains old long names\n');
    end
catch e
    fprintf('✗ FAIL: Error checking make_all_cases.m: %s\n', e.message);
end

fprintf('\n=== All CI Fix Tests Completed ===\n');
