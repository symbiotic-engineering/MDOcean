function [] = sweep_num_seeds(filename_uuid)
%SWEEP_NUM_SEEDS Sweep number of epsilon constraint seeds to see effect on
%runtime and hypervolume

    if nargin==0
        filename_uuid='';
    end

    num_seeds = 0:5:30;
    parfor i=1:length(num_seeds)
        [~,~,timeEpsConstraint(i),timeGradientFree(i),output] = pareto_search(filename_uuid,num_seeds(i))
        hypervolume(i) = output.volume;
        functioncountGradientFree(i) = output.funccount;
    end

    figure
    plot(num_seeds,timeEpsConstraint,...
         num_seeds,timeEpsConstraint+timeGradientFree,...
         num_seeds,hypervolume)
    legend('Time: Epsilon Constraint','Time: Total', 'Hypervolume')
end

