function intermed_result_struct = analysis_fcn(p,b)
    % Run Pareto search
    [r1_damping,r2_reactive] = damping_vs_reactive(p,b);

    intermed_result_struct.r1 = r1_damping;
    intermed_result_struct.r2 = r2_reactive;
end
