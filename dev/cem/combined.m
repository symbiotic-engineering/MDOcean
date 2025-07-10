clear
close all

p = parameters();
b = var_bounds();
X = [b.X_noms; 1];
[~, ~, ~, ~, val] = simulation(X,p);
w = val.w;
wave_amp = repmat(p.Hs,[1,size(p.JPD,2)]) / (2*sqrt(2));
F = val.gamma_f_over_rho_g * p.rho_w * p.g .* wave_amp;
RAO = false;       % true uses X/eta, false uses X/F_f
relative = false; % true uses X_u, false uses X_f
fudge = pi;       % phase shift correction

wn = 0.54;        % guess

phase_X = val.phase_X_f;
X       = val.X_f;

% get magnitude matrix
mag_matrix   = X ./ F;
angle_matrix = phase_X - val.gamma_phase_f - fudge;


fill_nan = false; 

% hole-filling (toggle with fill_nan)

if fill_nan
mag_matrix = fillmissing(mag_matrix, 'nearest');
angle_matrix = fillmissing(angle_matrix, 'nearest');
end

A      = pi/4 * b.D_f_nom^2 - b.D_s_nom^2;
K_h    = p.rho_w * p.g * A;
K_mult = 6;      % guess
K      = K_h * K_mult;
y_lab       = 'X/F';
y_lab_extra = 'K \omega_n^2';
norm        = K * wn^2;
use_db = true;

% Convert magnitude to dB, if flagged

mag_db = 20*log10(mag_matrix);

omega = w;


% define magnitude model (do not scale first if use_db is toggled off)

if use_db
mag_model = fittype(@(omega_n, k, zeta, w) ...
    20*log10(abs(((omega_n^2).*k) ./ (-w.^2 + 1i*2*zeta*omega_n.*w + omega_n^2))), ...
  'independent','w','coefficients',{'omega_n','k','zeta'});
else
mag_model = fittype(@(omega_n, k, zeta, w) ...
    abs(((omega_n^2).*k) ./ (-w.^2 + 1i*2*zeta*omega_n.*w + omega_n^2)), ...
  'independent','w','coefficients',{'omega_n','k','zeta'});   
end

%{
mag_model = fittype(@(omega_n, k, zeta, w) ...
    abs(((omega_n^2).*k) ./ (-w.^2 + 1i*2*zeta*omega_n.*w + omega_n^2)), ...
  'independent','w','coefficients',{'omega_n','k','zeta'});
%}


% define phase model

angle_model = fittype(@(omega_n, k, zeta, w) ...
    angle(((omega_n^2).*k) ./ (-w.^2 + 1i*2*zeta*omega_n.*w + omega_n^2))./pi, ...
  'independent','w','coefficients',{'omega_n','k','zeta'});

% preallocating states, since MATLAB was causing problems with differing
% array sizes after removing NaNs. should not throw warnings if fill_nan is
% toggled on

nStates = size(omega,1);
omega_n_fit       = nan(nStates,1);
zeta_fit          = nan(nStates,1);
k_fit             = nan(nStates,1);
omega_n_fit_phase = nan(nStates,1);
zeta_fit_phase    = nan(nStates,1);
k_fit_phase       = nan(nStates,1);

avgs = mean(omega,1,'omitnan');
lwb  = 7;  % high‐frequency cutoff index
tol  = 4;  % min points for fitting



plot_omegas = omega(~isnan(omega));
plot_omegas = unique(plot_omegas,'stable');

for i = 1:nStates

    if use_db
    mags   = mag_db(i,:);
    else
    mags = mag_matrix(i,:);
    end

    omegas = omega(i,:);
    angles = angle_matrix(i,:);
    

    col = [(i-1)/(nStates-1)  0  1-(i-1)/(nStates-1)];
    
    % PLOT MAGNITUDES
    nexttile(1)

    loglog(omegas/wn, mags, '*-', 'Color',col, 'DisplayName',sprintf('H_s=%.2f',p.Hs(i)))

    
    ylabel(['Magnitude |' y_lab '| ' y_lab_extra ' (-)'])
    hold on
    

    % PLOT PHASES
    nexttile(2)
    semilogx(omegas/wn, angles./pi, '*-', 'Color',col, 'LineWidth',2.0)
    xlabel('\omega/\omega_n (-)')
    ylabel(['Phase \angle(' y_lab ') / \pi'])
    hold on
    
    % remove high frequencies (lwb) and clear NaNs (unused if fill_nan is
    % true)
    mags_c   = mags(lwb:end);
    om_c     = omegas(lwb:end);
    ang_c    = angles(lwb:end);
    om_c(isnan(om_c)) = avgs(isnan(om_c));
    
    om_c_mag   = om_c(~isnan(mags_c));
    om_c_ang   = om_c(~isnan(ang_c));
    mags_c     = mags_c(~isnan(mags_c));
    ang_c      = ang_c(~isnan(ang_c));
    
    % fit magnitudes, get fit params
    if numel(mags_c) > tol

        [mf, gof_m] = fit(om_c_mag.', mags_c.', mag_model, 'StartPoint',[0.5,1e-6,0.1]);


        omega_n_fit(i) = mf.omega_n;
        zeta_fit(i)    = mf.zeta;
        k_fit(i)       = mf.k;
        
        fprintf('Magnitude fit H_s=%.2f: R^2 = %.4f\n', p.Hs(i), gof_m.rsquare);

        w_fit = logspace(log10(min(omega(:),[],'omitnan')), ...
                         log10(max(omega(:),[],'omitnan')), 500);
        mag_fit_vals = 20*log10(abs((mf.k.*mf.omega_n^2) ./ ...
                          (-w_fit.^2 + 1i*2*mf.zeta*mf.omega_n.*w_fit + mf.omega_n^2)));
        nexttile(1)
        plot(w_fit/wn, mag_fit_vals,'-','Color',col,'DisplayName',sprintf('H_s=%.2f',p.Hs(i)),'Color',[((i-1)/nStates) 0 (1 - (i-1)/nStates)], 'HandleVisibility','off','LineWidth',1.5)
    end
    
    % fit phasesm get fit params
    if numel(ang_c) > tol
        [pf, gof_p] = fit(om_c_ang.', ang_c.', angle_model, 'StartPoint',[0.5,1e-6,0.1]);
        omega_n_fit_phase(i) = pf.omega_n;
        zeta_fit_phase(i)    = pf.zeta;
        k_fit_phase(i)       = pf.k;
        
        fprintf('Phase fit H_s=%.2f: R^2 = %.4f\n', p.Hs(i), gof_p.rsquare);

        phase_fit_vals = angle((pf.k.*pf.omega_n^2) ./ ...
                          (-w_fit.^2 + 1i*2*pf.zeta*pf.omega_n.*w_fit + pf.omega_n^2)) ./ pi;
        nexttile(2)

        plot(w_fit/wn, phase_fit_vals,'-','Color',col,'DisplayName',sprintf('\\zeta=%.3f',pf.zeta),'Color',[((i-1)/nStates) 0 (1 - (i-1)/nStates)],'HandleVisibility','off','LineWidth',1.5)
    end
end



% plot formatting

%nexttile(1); legend('Location','eastoutside');


%custom color scheme
rng = linspace(0, (nStates-1)/(nStates-1), 15);
color_scheme = [rng; zeros(size(rng)); flip(rng)].';

nexttile(1);
cb1 = colorbar('eastoutside');
cb1.Label.String = 'H_s (m)';
colormap(color_scheme);
clim([ min(p.Hs)  max(p.Hs) ]); 
cb1.Label.FontSize  = 20;

%nexttile(2); legend('Location','eastoutside');

nexttile(2);
cb2 = colorbar('eastoutside');
colormap(color_scheme);
clim([ min(p.Hs)  max(p.Hs) ]); 
cb2.Label.String = 'H_s (m)';
cb2.Label.FontSize  = 20;

fitResults = table(p.Hs(:), omega_n_fit, zeta_fit, k_fit, ...
                   'VariableNames',{'Hs','omega_n','zeta','k'});

%dummy plot
nexttile(1);
hold on;
hData = plot(nan, nan, '*-', 'Color', [0 0 0], 'LineWidth', 1,   'HandleVisibility','on');
hFit  = plot(nan, nan, '-',  'Color', [0 0 0], 'LineWidth', 1.5, 'HandleVisibility','on');
hold off;
legend([hData, hFit], {'MDOcean','Fit'}, 'Location','southwest');

%{
nexttile(2);
hold on;
hData = plot(nan, nan, '*-', 'Color', [0 0 0], 'LineWidth', 1,   'HandleVisibility','on');
hFit  = plot(nan, nan, '-',  'Color', [0 0 0], 'LineWidth', 1.5, 'HandleVisibility','on');
hold off;
legend([hData, hFit], {'MDOcean','Fit'}, 'Location','southeast');
%}

% indicate area of unused high frequencies
nexttile(1)
xregion(plot_omegas(lwb)/wn, plot_omegas(1)/wn, "FaceColor","black","FaceAlpha",0.1,'HandleVisibility','off')


nexttile(2)
xregion(plot_omegas(lwb)/wn, plot_omegas(1)/wn, "FaceColor","black","FaceAlpha",0.1,'HandleVisibility','off')



disp(fitResults)

fitPhaseResults = table(p.Hs(:), omega_n_fit_phase, zeta_fit_phase, k_fit_phase, ...
                        'VariableNames',{'Hs','omega_n','zeta','k'});
disp(fitPhaseResults)

plot(NaN,NaN,'k*-','DisplayName','RM3')

improvePlot