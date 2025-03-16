%% numerically
function b_inf_numeric()
    a1_num = 6/100;
    a2_num = 20/100;
    d1_num = 75/100;
    d2_num = 60/100;
    h_num = 1;
    
    m0h_max_zero = find_m0h_max(a1_num,a2_num,d1_num,d2_num,h_num);
    
    plot_b_approx(a1_num,a2_num,d1_num,d2_num,h_num,m0h_max_zero)

end

function m0h_max_zero = find_m0h_max(a1_num,a2_num,d1_num,d2_num,h_num)
    % find m0h_max numerically
    d2s = linspace(.15,d1_num-eps,20);
    m0h_maxs_nan = acosh(realmax) ./ (1-d2s/h_num);
    m0h_max_zero = acosh(realmax)/2;
    m0h_max_asymp = acosh(realmax) ./ d2s/h_num;
    m0h_max_actual = zeros(length(d2s),1);
    m0h_max_actual_asymptotic = zeros(length(d2s),1);
    for j=1:length(d2s)
        upper = max(m0h_maxs_nan(j),m0h_max_asymp(j));
        lower = min(m0h_maxs_nan(j),m0h_max_asymp(j));
        m0 = linspace(.9*lower, upper*1.1, 800);

        b_entry = zeros(length(m0),1);
        b_high_freq = zeros(length(m0),1);
        for i=1:length(m0)
            [b_entry(i),b_high_freq(i)] = get_b(m0(i),h_num,a1_num,a2_num,d1_num,d2s(j));
        end
        m0h_max_actual(j) = max(m0(isfinite(b_entry)));
        m0h_max_actual_asymptotic(j) = max(m0(b_high_freq~=0));
    end
    
    figure
    yyaxis left
    plot(d2s,m0h_max_actual/acosh(realmax))
    hold on 
    plot(d2s,m0h_max_actual./m0h_maxs_nan')
    xlabel('d2/h')
    ylabel('normalized m_0h max for NaN')
    yyaxis right
    ylabel('normalized m_0h max for asymptotic 0')
    plot(d2s,m0h_max_actual_asymptotic/acosh(realmax))
    plot(d2s,m0h_max_actual_asymptotic./m0h_max_asymp')
    legend('normalize by acosh(realmax)','normalize by acosh(realmax)/(1-d_2/h)',...
        'normalize by acosh(realmax)','normalize by acosh(realmax)/(d_2/h)')
end

function plot_b_approx(a1_num,a2_num,d1_num,d2_num,h_num,m0h_max_zero)
    % check that my formula for b_high_freq is correct
    
    m0h_max_nan = acosh(realmax) / (1-d2_num/h_num);
%     lambert_arg = ( a2_num / (realmin * (1-d2_num/h_num)) )^2 * d2_num/h_num;
    m0h_max_asymptotic_zero = acosh(realmax) ./ d2_num/h_num;
    
    m0 = sort([logspace(log10(m0h_max_zero/100),log10(m0h_max_asymptotic_zero*2)) ...
                m0h_max_nan m0h_max_zero m0h_max_zero+eps ...
                m0h_max_asymptotic_zero m0h_max_asymptotic_zero+eps]);
    b_entry = zeros(length(m0),1);
    b_high_freq = zeros(length(m0),1);
    
    for i=1:length(m0)
        [b_entry(i),b_high_freq(i)] = get_b(m0(i),h_num,a1_num,a2_num,d1_num,d2_num);
    end
    
    display([m0' b_entry b_high_freq])
    
    b_both = [b_entry; b_high_freq];
    near_zero_b = max(nonzeros(b_both))/10;
    b_entry(b_entry==0) = near_zero_b; % so zero shows up on log-log plot
    b_high_freq(b_high_freq==0) = near_zero_b;
    
    figure
    loglog(m0/acosh(realmax),b_entry,'DisplayName','b')
    hold on
    loglog(m0/acosh(realmax),b_high_freq,'--','DisplayName','Asymptotic b')
    
    ylim([min(b_both) near_zero_b])
%     plot(m0h_max_nan/acosh(realmax) *[1 1],ylim,'k--','DisplayName','1/(1-d_2/h)')
    plot(m0h_max_zero/acosh(realmax)*[1 1],ylim,'k:','DisplayName','1/2')
    plot(m0h_max_asymptotic_zero/acosh(realmax)*[1 1],ylim,'k-.','DisplayName','1/(d_2/h)')
    
    legend
    xlabel('m_0h/acosh(realmax)')
    ylabel('b_{N+2M+1}')
    improvePlot
    
    % conclusion: it's not correct, but it goes to zero before it goes to nan,
    % so I'm ok to just zero it when above the m0h_max_nan threshold.
end

function [b_entry,b_high_freq] = get_b(m0_num,h_num,a1_num,a2_num,d1_num,d2_num)
    K_num = 10;
    m_k_cell = get_m_k(m0_num,h_num,K_num);
    
    [A_num, b_num, c_num, c_0_num] = A_b_c_matrix_N10_M10_K10_heaving_outer(a1_num, a2_num, d1_num, d2_num,...
                                h_num, m0_num, m_k_cell{:});
    b_entry = b_num(end-K_num);
    b_high_freq = -a2_num / (h_num - d2_num) * sqrt(h_num / (2 * m0_num)) * exp(-d2_num * m0_num);
end

function m_k_cell = get_m_k(m0_num,h_num,K_num)
% copy-pasted from run_MEEM/compute_eigen_hydro_coeffs
    m_k_num = zeros(1,K_num);
    % using tand instead of tan because finite precision of pi means
    % inconsistent behaviour around tan(pi/2)
    eqn = @(m_k_h_deg) (m_k_h_deg * pi/180) .* tand(m_k_h_deg) + m0_num * h_num * tanh(m0_num * h_num);

    for k_num = 1:K_num
        bounds = 180 * [k_num-1/2, k_num];
        % apply tweak determined experimentally to be the smallest tweak
        % that doesn't return +-Inf - see p17 of notebook
        expo = floor( log(bounds(1)) / log(2) );
        bound_tweak = eps * (2^(expo - 1) + 1);
        bounds(1) = bounds(1) + bound_tweak;

        m_k_h_deg = fzero(eqn, bounds);
        m_k_num(k_num) = m_k_h_deg * pi/180 / h_num;
    end    

    m_k_cell = num2cell(m_k_num);
end