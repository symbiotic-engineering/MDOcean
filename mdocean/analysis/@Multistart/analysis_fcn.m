function intermed_result_struct = analysis_fcn(p,b)
    % Run multistart optimization analysis
    [X_opt,objs,flags,x0s,num_runs] = gradient_mult_x0(p,b);
    
    % Store results for post-processing
    intermed_result_struct.p = p;
    intermed_result_struct.b = b;
    intermed_result_struct.X_opt = X_opt;
    intermed_result_struct.objs = objs;
    intermed_result_struct.flags = flags;
    intermed_result_struct.x0s = x0s;
    intermed_result_struct.num_runs = num_runs;
end