%{
lwb = 6;

omega = [1.142, 0.967, 0.838, 0.739, 0.661, 0.598, 0.546, ...
         0.503, 0.465, 0.433, 0.405, 0.381]; % Frequency (rad/s)
omega = omega(lwb:end);

mag_lin = [7.19e-07, 4.90e-07, 4.16e-07, 3.82e-07, 3.79e-07, ...
           4.51e-07, 8.72e-07, 1.00e-06, 6.68e-07, 5.36e-07, ...
           4.75e-07, 4.42e-07]; % Magnitude (linear)
mag_lin = mag_lin(lwb:end);

phase = [-0.009174, -0.041284, -0.077982, -0.114679, ...
         -0.169725, -0.220596, -0.270222, -0.310012, ...
         -0.348858, -0.387704, -0.428700, -0.471168]; % Phase (radians)
phase = phase(lwb:end);
%}

p = parameters();
b = var_bounds();
X = [b.X_noms; 1];
[~, ~, ~, ~, val] = simulation(X,p);
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

% setup fit model
mag_model = fittype(@(omega_n, k, zeta, w) ...
    20*log10(abs(((omega_n^2) .* k) ./ (-w.^2 + 1i*2*zeta*omega_n.*w + omega_n^2))), ...
    'independent', 'w', 'coefficients', {'omega_n', 'k', 'zeta'});

% determine cutoff for used data range, to avoid high frequency inaccuracy

%{
loglog(w(i,:)/wn, mag_norm, '*-', ...
        'DisplayName',['H_s=' num2str(p.Hs(i))],'Color',col)
    hold on
%}

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


for i = 1:size(omega,1)
    mags = mag_db(i, :)
    omegas = omega(i, :)

    
%{
% input data
[mag_fit, gof_mag] = fit(omegas.', mags.', mag_model, ...
    'StartPoint', [0.5, 1e-6, 0.1]);

w_fit = logspace(log10(min(omega)), log10(max(omega)), 500);
mag_fit_vals = 20*log10(abs((mag_fit.k .* mag_fit.omega_n^2) ./ ...
    (-w_fit.^2 + 1i*2*mag_fit.zeta*mag_fit.omega_n.*w_fit + mag_fit.omega_n^2)));
%}


    %{
    color_frac = (i-1)/(length(p.Hs)-1);
    col = [color_frac 0 1-color_frac];
    mag_norm = mag_data * norm;
    
    nexttile(1)
    loglog(w(i,:)/wn, mag_norm, '*-', ...
        'DisplayName',['H_s=' num2str(p.Hs(i))],'Color',col)
    hold on

    nexttile(2)
    h = semilogx(w(i,:)/wn, angle_data/pi,'*-','Color',col);
    h.Annotation.LegendInformation.IconDisplayStyle = 'off'; 
    hold on
    %}

    color_frac = (i-1)/(length(p.Hs)-1);
    col = [color_frac 0 1-color_frac];
    mag_norm = mag_data * norm;
    
    nexttile(1)
    loglog(omegas/wn, mags, '*-', ...
        'DisplayName',['H_s=' num2str(p.Hs(i))],'Color',col)
    hold on

    
    nexttile(2)
    h = semilogx(omegas, angle_data/pi,'*-','Color',col);
    h.Annotation.LegendInformation.IconDisplayStyle = 'off'; 
    hold on
    



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
[mag_fit, gof_mag] = fit(omega.', mag_db.', mag_model, ...
    'StartPoint', [0.5, 1e-6, 0.1]);
%{
% === FIT PHASE MODEL ===
phase_model = fittype(@(omega_n, k, zeta, w) ...
    angle((omega_n^2) ./ (-w.^2 + 1i*2*zeta*omega_n.*w + omega_n^2)), ...
    'independent', 'w', 'coefficients', {'omega_n', 'k', 'zeta'});

[phase_fit, gof_phase] = fit(omega.', phase.', phase_model, ...
    'StartPoint', [mag_fit.omega_n, mag_fit.zeta]);
%}
% === EVALUATE FITTED CURVES ===
w_fit = logspace(log10(min(omega)), log10(max(omega)), 500);
mag_fit_vals = 20*log10(abs((mag_fit.k .* mag_fit.omega_n^2) ./ ...
    (-w_fit.^2 + 1i*2*mag_fit.zeta*mag_fit.omega_n.*w_fit + mag_fit.omega_n^2)));
%{
phase_fit_vals = angle((phase_fit.omega_n^2) ./ ...
    (-w_fit.^2 + 1i*2*phase_fit.zeta*phase_fit.omega_n.*w_fit + phase_fit.omega_n^2));
%}
% === PLOTS ===
figure;
subplot(2,1,1);
semilogx(omega, mag_db, 'r.', w_fit, mag_fit_vals, 'b', 'LineWidth', 2);
xlabel('Frequency (rad/s)');
ylabel('Magnitude (dB)');
title('Magnitude Fit');
legend('Data', 'fit');
grid on;
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
fprintf('Fitted omega_n: %.4f rad/s\n', mag_fit.omega_n);
fprintf('Fitted k:    %.4f\n', mag_fit.k);
fprintf('Fitted zeta:    %.4f\n', mag_fit.zeta);


disp('Goodness of Fit (Magnitude):'); disp(gof_mag);
%{
disp('Goodness of Fit (Phase):'); disp(gof_phase);
%}