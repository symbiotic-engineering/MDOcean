function intermed_result_struct = analysis_fcn(p,b)
    % Run comparison analysis
    [Xs,vals] = compare_run(p,b);

    % Store results and returned figure handles for post-processing
    intermed_result_struct.p = p;
    intermed_result_struct.b = b;
    intermed_result_struct.Xs = Xs;
    intermed_result_struct.vals = vals;
end
