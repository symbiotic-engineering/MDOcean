clear
close all

durations = 1;
DP = 1;
E = 1:10;

for i_tau = 1:length(durations)
    tau = durations(i_tau);
    for i_DP = 1:length(DP)
        JPD = 1;
        P_wec_DP = 1;
        P_wec_tau_DP = 1;
        for i = 1:length(E)
            P_lim = E / tau;
            P_load = 1;
            p_empty = 1;
            p_full = 1;
            P_load = P_load * 1;
            cv = 1;
        end
        figure
        plot(E,cv)
        DP_star = DP(i_DP) : DP(end);
        for i_DPstar = 1:length(DP_star)
            E_DPstar = 1;
            c = 1;
        end
    end
    c_min = 1;
end

% get power PDF for min LCOE design

% plot pareto front
pareto_curve_heuristics()
figure(2)

