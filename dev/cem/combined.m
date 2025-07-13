clear
close all

%% choose settings
RAO = false;       % true uses X/eta, false uses X/F_f
relative = false; % true uses X_u, false uses X_f
phase_uses_mag_fit = false;
fill_nan = false; 
use_db_for_fit = true;
use_db_for_plot = false;

%% define data
p = parameters();
p.JPD(p.JPD==0) = 1;

b = var_bounds();
X = [b.X_noms; 1];
[~, ~, ~, ~, val] = simulation(X,p);
w = val.w;
wave_amp = repmat(p.Hs,[1,size(p.JPD,2)]) / (2*sqrt(2));
F = val.gamma_f_over_rho_g * p.rho_w * p.g .* wave_amp;

if phase_uses_mag_fit
    fudge = 0;       % phase shift correction
else
    fudge = pi;
end

% wn = 0.54;        % guess

phase_X = val.phase_X_f;
X       = val.X_f;

% get magnitude matrix

mag_matrix   = X ./ F;
angle_matrix = phase_X - val.gamma_phase_f - fudge;

% hole-filling (toggle with fill_nan)

if fill_nan
    mag_matrix = fillmissing(mag_matrix, 'nearest');
    angle_matrix = fillmissing(angle_matrix, 'nearest');
end

% A      = pi/4 * (b.D_f_nom^2 - b.D_s_nom^2);
% K_h    = p.rho_w * p.g * A;
% K_mult = 6;      % guess
% K      = K_h * K_mult;
y_lab       = 'X/F';
omega = w;

%% second order model definition
sec_order_fn = @(w_n, k, zeta, w) 1/k * 1./(1 - (w/w_n).^2 + 1i * 2*zeta*w/w_n);
mag_fn       = @(w_n, k, zeta, w) abs( sec_order_fn(w_n, k, zeta, w) );
mag_db_fn    = @(w_n, k, zeta, w) 20*log10( mag_fn(w_n, k, zeta, w) );
phase_fn     = @(w_n, k, zeta, w) 1/pi * angle( sec_order_fn(w_n, k, zeta, w) );

%% for magnitude fit and plot, convert magnitude to dB, if flagged
mag_db = 20*log10(mag_matrix);
scale_mag_data_for_fit = 1e6; % to keep values O(1) numerically (helps convergence)
if use_db_for_fit
    mag_data_for_fit = mag_db + 20*log10(scale_mag_data_for_fit);
    mag_fn_for_fit = mag_db_fn;
else
    mag_data_for_fit = mag_matrix * scale_mag_data_for_fit;
    mag_fn_for_fit = mag_fn;
end

if use_db_for_plot
    mag_data_for_plot = mag_db;
    mag_fn_for_plot = mag_db_fn;
    mag_plot_fn = @(varargin)semilogx(varargin{:}); % y axis already log so don't take log again
    y_add = ', dB';
else
    mag_data_for_plot = mag_matrix;
    mag_fn_for_plot = mag_fn;
    mag_plot_fn = @(varargin)loglog(varargin{:});
    y_add = '';
end

%% define mag and phase fit models
mag_model = fittype(mag_fn_for_fit, ...
    'independent','w','coefficients',{'w_n','k','zeta'});

angle_model = fittype(phase_fn, ...
  'independent','w','coefficients',{'w_n','k','zeta'});

%% preallocating states
% since MATLAB was causing problems with differing
% array sizes after removing NaNs. should not throw warnings if fill_nan is
% toggled on

nWaveHeights = size(omega,1);
omega_n_fit       = nan(nWaveHeights,1);
zeta_fit          = nan(nWaveHeights,1);
k_fit             = nan(nWaveHeights,1);
omega_n_fit_phase = nan(nWaveHeights,1);
zeta_fit_phase    = nan(nWaveHeights,1);
k_fit_phase       = nan(nWaveHeights,1);

avgs = mean(omega,1,'omitnan');
lwb  = 8;  % high‐frequency cutoff index
tol  = 4;  % min points for fitting

%% loop through wave heights
for i = 1:nWaveHeights
    omegas = omega(i,:);
    angles = angle_matrix(i,:)/pi;
    
    red = (i-1)/(nWaveHeights-1);
    green = 0;
    blue = 1 - red;
    col = [red green blue];
    
    % PLOT MAGNITUDE DATA
    nexttile(1)
    mag_plot_fn(omegas, mag_data_for_plot(i,:), '*--', 'Color',col, 'DisplayName',sprintf('H_s=%.2f',p.Hs(i)))
    ylabel(['Magnitude |' y_lab '|' y_add])
    hold on

    % PLOT PHASE DATA
    nexttile(2)
    semilogx(omegas, angles, '*--', 'Color',col,'LineWidth',2.0)
    xlabel('Frequency \omega (rad/s)')
    ylabel(['Phase \angle(' y_lab ') / \pi'])
    hold on
    
    % remove high frequencies (lwb)
    mags_c   = mag_data_for_fit(i,lwb:end);
    om_c     = omegas(lwb:end);
    ang_c    = angles(lwb:end);
    om_c(isnan(om_c)) = avgs(isnan(om_c));
    
    % clear NaNs (unused if fill_nan is true)
    idx_nan = isnan(mags_c);
    om_c   = om_c(~idx_nan);
    mags_c = mags_c(~idx_nan);
    ang_c  = ang_c(~idx_nan);
    
    % fit magnitudes, get fit params
    if numel(mags_c) > tol

        [mf, gof_m] = fit(om_c.', mags_c.', mag_model, 'StartPoint',[0.5,5,0.05]);

        omega_n_fit(i) = mf.w_n;
        zeta_fit(i)    = mf.zeta;
        k_fit(i)       = mf.k * scale_mag_data_for_fit;
        
        fprintf('Magnitude fit H_s=%.2f: R^2 = %.4f\n', p.Hs(i), gof_m.rsquare);

        w_fit = logspace(log10(min(omega(:),[],'omitnan')), ...
                         log10(max(omega(:),[],'omitnan')), 500);
        mag_fit_vals = mag_fn_for_plot(mf.w_n, mf.k * scale_mag_data_for_fit, mf.zeta, w_fit);

        nexttile(1)
        mag_plot_fn(w_fit, mag_fit_vals,'-','Color',[col .5],...
            'DisplayName',sprintf('H_s=%.2f',p.Hs(i)), ...
            'HandleVisibility','off','LineWidth',1.5)

        ymin = 1e-7;
        if use_db_for_plot
            ymin = 20*log10(ymin);
        end
        ylim([ymin max(mag_data_for_plot,[],'all')])

    end
    

    % fit phases, get fit params
    if numel(ang_c) > tol
        if phase_uses_mag_fit
            pf = mf;
            gof_p = gof_m;
            omega_n_fit_phase(i) = NaN;
            zeta_fit_phase(i)    = NaN;
            k_fit_phase(i)       = NaN;
        else
            [pf, gof_p] = fit(om_c.', ang_c.', angle_model, 'StartPoint',[0.5,5,0.05]);
            omega_n_fit_phase(i) = pf.w_n;
            zeta_fit_phase(i)    = pf.zeta;
            k_fit_phase(i)       = pf.k * scale_mag_data_for_fit;

        end
        
        fprintf('Phase fit H_s=%.2f: R^2 = %.4f\n', p.Hs(i), gof_p.rsquare);

        phase_fit_vals = angle_model(pf.w_n, pf.k, pf.zeta, w_fit);
        nexttile(2)

        semilogx(w_fit, phase_fit_vals,'-','Color',[col 0.5],...
            'DisplayName',sprintf('\\zeta=%.3f',pf.zeta),...
            'HandleVisibility','off','LineWidth',1.5)
    end
end



%% plot formatting

%nexttile(1); legend('Location','eastoutside');


%custom color scheme
rng = linspace(0, (nWaveHeights-1)/(nWaveHeights-1), 15);
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
%xregion(plot_omegas(lwb), plot_omegas(1)/wn, "FaceColor","black","FaceAlpha",0.1,'HandleVisibility','off')


nexttile(2)
%xregion(plot_omegas(lwb), plot_omegas(1)/wn, "FaceColor","black","FaceAlpha",0.1,'HandleVisibility','off')

plot(NaN,NaN,'k*--','DisplayName','RM3')

improvePlot
set(gcf(),"Position",[100 100 600 680])

%% display of fit tables

fitResults = table(p.Hs(:), omega_n_fit, zeta_fit, k_fit, ...
                   'VariableNames',{'Hs','omega_n','zeta','k'});

disp(fitResults)

fitPhaseResults = table(p.Hs(:), omega_n_fit_phase, zeta_fit_phase, k_fit_phase, ...
                        'VariableNames',{'Hs','omega_n','zeta','k'});
disp(fitPhaseResults)

%% Dependence of mag and phase fits on Hs
figure
t = tiledlayout(3,1);

nexttile
plot(p.Hs(:), omega_n_fit, 'g')
hold on
plot(p.Hs(:), omega_n_fit_phase,'m-.')
ylabel({'Natural','Frequency','$\omega_n$'},'Interpreter','latex')
grid on 

nexttile
plot(p.Hs(:), zeta_fit, 'g')
hold on
plot(p.Hs(:), zeta_fit_phase,'m-.')
ylabel({'Damping','Ratio','$\zeta$'},'Interpreter','latex')
grid on 

nexttile
h = plot(p.Hs(:), k_fit, 'g');
h.Annotation.LegendInformation.IconDisplayStyle = 'off'; % no legend here
ylabel('Stiffness $K$','Interpreter','latex')
xlabel('Wave Height $H_s$ (m)','Interpreter','latex')
hold on
%plot(p.Hs(:), k_fit_phase,'m-.') K fit for phase is meaningless because you
%can't get any information about K from phase, only from magnitude
grid on 

plot(NaN,NaN,'g','DisplayName','Magnitude Fit')
plot(NaN,NaN,'m-.','DisplayName','Phase Fit')
legend
improvePlot

