%% symbolically
% from A_b_c_matrix_N10_M10_K10_heaving_outer, tracing back the top (k=0)
% term of the bottom block of the b-vector
clear
syms a1 a2 d1 d2 h m0 real positive

t16 = d2.*2.0;
t19 = h.*2.0;
t35 = m0.*t19;
t51 = -h;
t52 = -t19;
t53 = 1.0./h;
t54 = 1.0./m0;
t81 = sinh(t35);
t161 = d2+t51;
t244 = m0.*t161;
t255 = t16+t52;
t261 = sinh(t244);
t295 = 1.0./t255;

b0 = -a2.*t54.*t261.*t295.*1.0./sqrt((t53.*t54.*t81)./4.0+1.0./2.0);
b0 = subs(expand(b0),sinh(m0*h), cosh(m0*h)*tanh(m0*h));
[num,den] = numden(simplifyFraction(b0));
new_num = simplifyFraction(num/cosh(h*m0));
new_den = sqrt(partfrac(den^2/cosh(h*m0)^2));
pretty(new_num/new_den)
b0_approx = subs(new_num/new_den, {tanh(h*m0),1/cosh(h*m0)^2},{1,0});
pretty(simplify(b0_approx))

%% numerically

close all
a1_num = 6/100;
a2_num = 20/100;
d1_num = 75/100;
d2_num = 60/100;
h_num = 1;

% find m0h_max numerically
b_entry = zeros(length(m0),1);
b_high_freq = zeros(length(m0),1);
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

%% check that my formula for b_high_freq is correct

m0h_max_nan = acosh(realmax) / (1-d2_num/h_num);
lambert_arg = ( a2_num / (realmin * (1-d2_num/h_num)) )^2 * d2_num/h_num;
m0h_max_asymptotic_zero = acosh(realmax) ./ d2_num/h_num;

m0 = sort([logspace(log10(m0h_max_zero/500),log10(m0h_max_asymptotic_zero*2)) ...
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
loglog(m0/acosh(realmax),b_entry,m0/acosh(realmax),b_high_freq,'--')
hold on

ylim([min(b_both) near_zero_b])
plot(m0h_max_nan/acosh(realmax) *[1 1],ylim,'k--')
plot(m0h_max_zero/acosh(realmax)*[1 1],ylim,'k:')
plot(m0h_max_asymptotic_zero/acosh(realmax)*[1 1],ylim,'k-.')

legend('b','Asymptotic b','1/(1-d_2/h)','1/2','1/(d_2/h)')
xlabel('m_0h/acosh(realmax)')
ylabel('b_{34}')
improvePlot

% conclusion: it's not correct, but it goes to zero before it goes to nan,
% so I'm ok to just zero it when above the m0h_max_nan threshold.

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