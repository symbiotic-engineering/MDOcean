clear
close all

p = parameters();
b = var_bounds();
X = [b.X_noms; 1];
[~, ~, ~, val] = simulation(X,p);
w = val.w;
wave_amp = repmat(p.Hs,[1,size(p.JPD,2)]) / (2*sqrt(2));
F = val.gamma_f_over_rho_g * p.rho_w * p.g .* wave_amp;
RAO = true; % true uses X/eta, false uses X/F_f
relative = false; % true uses X_u, false uses X_f
fudge = pi; % unsure why this is necessary, otherwise I get a -180deg phase shift as if there were a negative sign

wn = 0.54 % guess

phase_X = val.phase_X_f;
X = val.X_f;

% normalize - no RAO
    mag_matrix = X ./ F;
    angle_matrix = phase_X - val.gamma_phase_f - fudge;
    %A = pi/4 * b.D_f_nom^2 - (p.D_f_in_over_D_s * b.D_s_nom)^2;
    A = pi/4 * b.D_f_nom^2 - b.D_s_nom^2;
    K_h = p.rho_w * p.g * A;
    K_mult = 6; % guess for now
    K = K_h * K_mult; % no pto stiffness
    y_lab = 'X/F';
    y_lab_extra = 'K\omega_n^2';
    norm = K * wn^2;


% Convert magnitude to dB
mag_db = 20 * log10(mag_matrix); 

% pull data from val/rename
omega = w;

% setup fit model for magnitude
mag_model = fittype(@(omega_n, k, zeta, w) ...
    20*log10(abs(((omega_n^2) .* k) ./ (-w.^2 + 1i*2*zeta*omega_n.*w + omega_n^2))), ...
    'independent', 'w', 'coefficients', {'omega_n', 'k', 'zeta'});

% setup fit model for angle
angle_model = fittype(@(omega_n, k, zeta, w) ...
    20*log10(abs(((omega_n^2) .* k) ./ (-w.^2 + 1i*2*zeta*omega_n.*w + omega_n^2))), ...
    'independent', 'w', 'coefficients', {'omega_n', 'k', 'zeta'});



%{
mag_data   = mag_matrix(5,:);
angle_data = angle_matrix(5,:);
TF_data = mag_data .* exp(1i * angle_data);
TF_w = w(5,:);

mag_nan   = mean(mag_matrix(:,isnan(TF_data)),'omitnan');
angle_nan = mean(angle_matrix(:,isnan(TF_data)),'omitnan');
TF_data(isnan(TF_data)) = mag_nan .* exp(1i * angle_nan);
unique_w = unique(w(:,isnan(TF_w)),'stable');
TF_w(isnan(TF_w)) = unique_w(~isnan(unique_w));

TF_frd = frd(TF_data(~isnan(TF_data)), TF_w(~isnan(TF_w)));
%}

% pass through and take averages for omega
avgs = mean(omega, 1, 'omitnan')


% cutoff for limiting points/minimum fit points required. This is due to a
% minimum of 3 points required for fittype to work, and since the highest
% <lwb> frequencies are discarded due to incorrect model assumptions, some
% of the fits have very little points to work off of which results in
% inaccurate models, the amount of points which can be controlled by <tol>
lwb = 5;
tol = 4;

for i = 1:size(omega,1)
    mags = mag_db(i, :);
    omegas = omega(i, :);
    angles = angle_matrix(i, :);


% plotting the data as connected points (solid indicates points used in
% fit, crosses indicate ignored)
    color_frac = (i-1)/(length(p.Hs)-1);
    col = [color_frac 0 1-color_frac];
    %mag_norm = mag_data * norm;
    
    nexttile(1)
    
    ylabel(['Magnitude |' y_lab '| ' y_lab_extra ' (-)'])
    %xlim([.5 omega_over_omega_n_max])

    wn = 0.54 % placeholder will be replaced later
    loglog(omegas, mags, '*-', ...
        'DisplayName',['H_s=' num2str(p.Hs(i))],'Color',col) % TODO: normalize by wn later
    hold on

    
    nexttile(2)

    xlabel('\omega/\omega_n (-)')
    ylabel(['Phase \angle(' y_lab ') / \pi'])

    h = semilogx(omegas, angles/pi,'*-','Color',col);
    h.Annotation.LegendInformation.IconDisplayStyle = 'off'; 
    hold on
    


    % remove NaN values + high freq values editable with lwb

    mags_clean = mags(1, lwb:end);
    omegas_clean = omegas(1, lwb:end);
    angles_clean = omegas(1, lwb:end);
    
    omegas_clean(isnan(omegas_clean)) = avgs(isnan(omegas_clean));
    omegas_clean_mag = omegas_clean(~isnan(mags_clean));
    omegas_clean_angle = omegas_clean(~isnan(angles_clean));

    mags_clean = mags_clean(~isnan(mags_clean));
    angles_clean = angles_clean(~isnan(angles_clean));

    

    % input fits and create fit models

    if size(mags_clean, 2) > tol

    [mag_fit, gof_mag] = fit(omegas_clean_mag.', mags_clean.', mag_model, ...
    'StartPoint', [0.5, 1e-6, 0.1]);

    % get results from fit
    
    w_fit = logspace(log10(min(omega, [], "all", "omitnan")), log10(max(omega,[], "all", "omitnan")), 500);

    mag_fit_vals = 20*log10(abs((mag_fit.k .* mag_fit.omega_n^2) ./ ...
    (-w_fit.^2 + 1i*2*mag_fit.zeta*mag_fit.omega_n.*w_fit + mag_fit.omega_n^2)));

    omega_over_omega_n_max = max(omega,[],'all'); % TODO: implement normalization by wn later

    %fplot(omega,mag_fit_vals,[0 omega_over_omega_n_max],...
    %    'DisplayName',['\zeta=' num2str(zeta(i))],'Color',(i-1)/length(zeta)*[1 1 1])
    nexttile(1)
        plot(w_fit.',mag_fit_vals.',...
        'DisplayName',['\zeta=' num2str(i)],'Color',(i-1)/i*[1 1 1])
    
    else

    disp('failed - too few points')
    
    end
    
end
improvePlot
%{
figure;
subplot(2,1,1);
semilogx(omega, mag_db, 'r.', w_fit, mag_fit_vals, 'b', 'LineWidth', 2);

xlabel('Frequency (rad/s)');
ylabel('Magnitude (dB)');
title('Magnitude Fit');
legend('Data', 'fitx');
grid on;
%}




% input data
%[mag_fit, gof_mag] = fit(omega.', mag_db.', mag_model, ...
%    'StartPoint', [0.5, 1e-6, 0.1]);
%{
% === FIT PHASE MODEL ===
phase_model = fittype(@(omega_n, k, zeta, w) ...
    angle((omega_n^2) ./ (-w.^2 + 1i*2*zeta*omega_n.*w + omega_n^2)), ...
    'independent', 'w', 'coefficients', {'omega_n', 'k', 'zeta'});

[phase_fit, gof_phase] = fit(omega.', phase.', phase_model, ...
    'StartPoint', [mag_fit.omega_n, mag_fit.zeta]);
%}
% === EVALUATE FITTED CURVES ===
%w_fit = logspace(log10(min(omega)), log10(max(omega)), 500);
%mag_fit_vals = 20*log10(abs((mag_fit.k .* mag_fit.omega_n^2) ./ ...
%    (-w_fit.^2 + 1i*2*mag_fit.zeta*mag_fit.omega_n.*w_fit + mag_fit.omega_n^2)));
%{
phase_fit_vals = angle((phase_fit.omega_n^2) ./ ...
    (-w_fit.^2 + 1i*2*phase_fit.zeta*phase_fit.omega_n.*w_fit + phase_fit.omega_n^2));
%}
% === PLOTS ===
%figure;
%subplot(2,1,1);
%semilogx(omega, mag_db, 'r.', w_fit, mag_fit_vals, 'b', 'LineWidth', 2);
%xlabel('Frequency (rad/s)');
%ylabel('Magnitude (dB)');
%title('Magnitude Fit');
%legend('Data', 'fit');
%grid on;
%{
subplot(2,1,2);
semilogx(omega, phase, 'r.', w_fit, phase_fit_vals, 'b', 'LineWidth', 2);
xlabel('Frequency (rad/s)');
ylabel('Phase (rad)');
title('Phase Fit');
legend('Data', 'fit');
grid on;

%}

% === DISPLAY RESULTS ===
%fprintf('Fitted omega_n: %.4f rad/s\n', mag_fit.omega_n);
%fprintf('Fitted k:    %.4f\n', mag_fit.k);
%fprintf('Fitted zeta:    %.4f\n', mag_fit.zeta);


%disp('Goodness of Fit (Magnitude):'); disp(gof_mag);
%{
disp('Goodness of Fit (Phase):'); disp(gof_phase);
%}