% load data from MEEM (hardcoded, will replace later)
%dat = readmatrix('mdocean/mdocean/inputs/validation/MEEM_validation/damping.csv');
%w = dat(:,1);
%mag = dat(:,2);
%idx = w > 0;
%w = w(idx);
%mag = mag(idx);

p = parameters();
b = var_bounds();
X = [b.X_noms; 1];
[~, ~, ~, ~, val] = simulation(X,p);
w = val.w;
F = val.gamma_f_over_rho_g * p.rho_w * p.g * p.Hs / (2*sqrt(2)) ;
X = val.X_f;
mag_matrix = X ./ F;
mag = mag_matrix(1,:); % all freqs at lowest wave height



% create fit model
fit_model = fittype(@(zeta,wn,K,w) ...
          K.*wn.^2 ./ sqrt((wn.^2 - w.^2).^2 + (2*zeta*wn.*w).^2), ...
          'independent','w', ...
          'coefficients',{'zeta','wn','K'});
opt            = fitoptions(fit_model);
opt.StartPoint = [0.5  1.0   1.0];   % zeta, omega_n, gain K

% create fit
[fit_result, stats] = fit(w, mag, fit_model, opt);

% get parameters from fit
zeta = fit_result.zeta;
wn = fit_result.wn;
K = fit_result.K;

fprintf('zeta (damping ratio)  = %.4g\n', zeta);
fprintf('omega_n (rad/s)  = %.4g\n', wn);
fprintf('low frequency gain (K)  = %.4g', K);
fprintf('R-squared value = %.5f\n', stats.rsquare);

% fit curve
wfit = logspace(log10(min(w)), log10(max(w)), 400);
magfit = feval(fit_result, wfit);

% plotting
figure('Color','w');
semilogx(w,   20*log10(mag), '.', 'MarkerSize',5, 'DisplayName','Data'); hold on
semilogx(wfit,20*log10(magfit), 'LineWidth',1.4, 'DisplayName','2nd-order fit');
xlabel('Omega (rad/s)');
ylabel('Magnitude (dB)');
title('Bode Magnitude');
legend