function intermed_result_struct = analysis_fcn(~,~)
    % Run describing function demo and capture returned figure handles
    figs = sin_desc_fcn_demo();
    intermed_result_struct.created_figs = figs;
end
