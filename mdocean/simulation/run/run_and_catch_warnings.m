function [outputs, warning_hit, captured_text] = run_and_catch_warnings(fcn, warning_ids, nout)
% RUN_AND_CATCH_WARNINGS  Call a function once, silently capturing warnings.
%
% This function saves the current warning state (restoring it via onCleanup),
% configures warning display so that each warning's printed text contains its
% identifier string, then runs FCN exactly once inside evalc so that all
% Command-Window output — including warning messages — is captured rather than
% printed to the terminal.  After the call it checks the captured text for
% each identifier in WARNING_IDS.
%
% Warning states do NOT need to be in any particular state before this call;
% this function handles all promotion and restoration internally.
%
% Implementation detail: MATLAB's verbose warning mode adds the line
%   (Type "warning off <identifier>" ...)
% to every warning printout.  This makes the identifier searchable in the
% evalc-captured text without requiring warnings to be elevated to errors.
%
% :param fcn:         Zero-argument function handle, e.g. @() myfun(a, b).
% :param warning_ids: 1×M cell array of warning identifier strings to detect.
% :param nout:        Number of output arguments to request from fcn.
% :returns:
%   outputs       – 1×nout cell array of outputs from fcn.
%   warning_hit   – 1×M logical: true where that warning fired at least once.
%   captured_text – char vector of all text captured by evalc.

    % Save and restore ALL warning states on exit (regardless of errors).
    prev_states = warning;
    cleanup = onCleanup(@() warning(prev_states)); %#ok<NASGU>

    % Suppress backtrace (noisy stack-trace lines) and enable verbose mode.
    % Verbose mode appends "(Type "warning off <id>" ...)" to each warning,
    % making the identifier searchable in the evalc-captured output.
    warning('off', 'backtrace');
    warning('on',  'verbose');

    % Ensure every managed warning is enabled so it fires and is captured.
    % (It must not be 'off' or 'error'; we want it to print so evalc sees it.)
    for k = 1:numel(warning_ids)
        warning('on', warning_ids{k});
    end

    % Run the function exactly once.  evalc captures all Command-Window output,
    % including warning messages, so nothing reaches the terminal.
    outputs = cell(1, nout);
    captured_text = evalc('[outputs{1:nout}] = fcn();');

    % Identify which managed warning IDs appeared in the captured output.
    % Match against the verbose-mode hint text "warning off <id>" so we
    % detect only actual triggered warnings, not unrelated occurrences of
    % the identifier string in other output.
    warning_hit = false(1, numel(warning_ids));
    for k = 1:numel(warning_ids)
        pat = ['warning off ' regexptranslate('literalstr', warning_ids{k})];
        warning_hit(k) = ~isempty(regexp(captured_text, pat, 'once'));
    end
end
