function intermed_result_struct = analysis_fcn(~,~)
    % Run damping plate structural analysis
    figs = viz_damping_plate();
    
    intermed_result_struct.figs = [figs(end), figs(end-1), figs(end-7)];
end
