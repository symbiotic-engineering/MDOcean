function [outputs, warning_hit, caught_MEs] = run_and_catch_warnings(fcn, warning_ids, nout)
% RUN_AND_CATCH_WARNINGS  Run a function, silently capturing managed warnings.
%
% The warnings listed in WARNING_IDS must already be promoted to errors
% (via warning('error', id)) before calling this function.  Each time one
% of those warnings fires (as an error), it is recorded and then silenced
% (warning('off', id)) so that the call can be retried.  This repeats until
% the call completes without throwing a managed-warning error.
%
% :param fcn:         Zero-argument function handle, e.g. @() myfun(a, b).
% :param warning_ids: 1×M cell array of warning identifier strings to catch.
% :param nout:        Number of output arguments to request from fcn.
% :returns:
%   outputs     – 1×nout cell array of all outputs from fcn.
%   warning_hit – 1×M logical: true where that warning fired at least once.
%   caught_MEs  – cell array of every caught MException (may contain
%                 duplicates if the same warning fired on multiple retries).

    warning_hit = false(1, numel(warning_ids));
    caught_MEs  = {};
    while true
        try
            outputs = cell(1, nout);
            [outputs{1:nout}] = fcn();
            break
        catch ME
            warning_idx = find(strcmp(ME.identifier, warning_ids), 1);
            if isempty(warning_idx)
                rethrow(ME)
            end
            warning_hit(warning_idx) = true;
            caught_MEs{end+1} = ME; %#ok<AGROW>
            warning('off', ME.identifier)
        end
    end
end
