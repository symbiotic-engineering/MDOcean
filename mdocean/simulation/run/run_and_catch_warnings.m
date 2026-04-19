function [outputs, warning_hit, captured_text] = run_and_catch_warnings(fcn, warning_ids, nout)
% RUN_AND_CATCH_WARNINGS  Call a function once, silently capturing warnings.
%
% This function saves the current warning state (restoring it via onCleanup),
% configures warning display so that each warning's printed text contains its
% identifier string, then runs FCN exactly once inside evalc so that all
% Command-Window output — including warning messages — is captured rather than
% printed to the terminal.  After the call it counts occurrences of each
% identifier in WARNING_IDS in the captured text.
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
% :param warning_ids: 1×M cell array of warning identifier strings to detect,
%                     OR the string 'all' to detect all warnings of any
%                     identifier that fire during fcn.
% :param nout:        Number of output arguments to request from fcn.
% :returns:
%   outputs       – 1×nout cell array of outputs from fcn.
%   warning_hit   – When warning_ids is a cell: 1×M integer array counting
%                   how many times each warning fired.
%                   When warning_ids is 'all': containers.Map from warning
%                   identifier string to integer count (only IDs that fired
%                   are present as keys).
%   captured_text – char vector of all text captured by evalc.

    all_mode = ischar(warning_ids) && strcmp(warning_ids, 'all');

    % Save and restore ALL warning states on exit (regardless of errors).
    prev_states = warning;
    cleanup = onCleanup(@() warning(prev_states)); %#ok<NASGU>

    % Suppress backtrace (noisy stack-trace lines) and enable verbose mode.
    % Verbose mode appends "(Type "warning off <id>" ...)" to each warning,
    % making the identifier searchable in the evalc-captured output.
    warning('off', 'backtrace');
    warning('on',  'verbose');

    if all_mode
        % Enable all warnings so any that fire are captured.
        warning('on', 'all');
    else
        % Ensure every managed warning is enabled so it fires and is captured.
        % (It must not be 'off'; we want it to print so evalc sees it.)
        for k = 1:numel(warning_ids)
            warning('on', warning_ids{k});
        end
    end

    % Run the function exactly once.  evalc captures all Command-Window output,
    % including warning messages, so nothing reaches the terminal.
    outputs = cell(1, nout);
    captured_text = evalc('[outputs{1:nout}] = fcn();');

    if all_mode
        % Extract all unique warning identifiers and their occurrence counts.
        % The verbose pattern is: (Type "warning off <identifier>" ...)
        % Identifiers consist of word characters and colons (e.g. A:B:C).
        % The trailing " anchors the match to the end of the identifier in
        % MATLAB's verbose hint, preventing false matches from function output.
        toks = regexp(captured_text, 'warning off ([\w:]+)"', 'tokens');
        warning_hit = containers.Map('KeyType', 'char', 'ValueType', 'double');
        for m = 1:numel(toks)
            id = toks{m}{1};
            if isKey(warning_hit, id)
                warning_hit(id) = warning_hit(id) + 1;
            else
                warning_hit(id) = 1;
            end
        end
    else
        % Count occurrences of each managed warning ID in the captured output.
        % Match against the verbose-mode hint text 'warning off <id>"' so we
        % detect only actual triggered warnings, not unrelated string matches.
        % The trailing " anchors the match to MATLAB's exact verbose format.
        warning_hit = zeros(1, numel(warning_ids));
        for k = 1:numel(warning_ids)
            pat = ['warning off ' regexptranslate('literalstr', warning_ids{k}) '"'];
            warning_hit(k) = numel(regexp(captured_text, pat));
        end
    end
end
