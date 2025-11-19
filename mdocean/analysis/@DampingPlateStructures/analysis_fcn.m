function intermed_result_struct = analysis_fcn(~,~)
    % Run damping plate structural analysis
    addpath('../dev/structures/damping-plate');
    BoedoPrantilAnnularPlate()
    
    % Store figure numbers for post-processing
    n = gcf().Number;
    intermed_result_struct.figs = [figure(n), figure(n-1), figure(n-7)];
end
