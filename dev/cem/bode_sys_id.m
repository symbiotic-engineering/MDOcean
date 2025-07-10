clear
close all
f = figure;
tiledlayout(2,1)

% nominal design
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
wn = .54; % guess for now
if relative
    phase_X = val.phase_X_u;
    X = val.X_u;
else
    phase_X = val.phase_X_f;
    X = val.X_f;
end
if RAO
    mag_matrix = X ./ wave_amp;
    angle_matrix = phase_X - fudge;
    y_lab = 'X/\eta';
    y_lab_extra = '';
    norm = 1;
else
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
end


for i=1:length(p.Hs)
    mag_data = mag_matrix(i,:);
    angle_data = angle_matrix(i,:);
    
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

end

%% symbolic: standard second order system
zeta = sort([.03 .04 .05 .07  .15 .1 : .1 : .3]);
omega_n = 1 : .5 : 2;

%[ZETA,OMEGA_N] = meshgrid(zeta,omega_n);

syms K ZETA OMEGA_N omega_over_omega_n real positive
omega = omega_over_omega_n * OMEGA_N;
s = 1i * omega;
X_over_F = 1/K / (s^2 + 2.*ZETA.*OMEGA_N*s + OMEGA_N.^2);
X_over_F_norm = simplify(X_over_F * K * OMEGA_N.^2);
pretty(X_over_F_norm)

mag = abs(X_over_F_norm);
phase = angle(X_over_F_norm);

omega_over_omega_n_max = max(w/wn,[],'all');
for i=1:length(zeta)
    nexttile(1)
    fplot(omega_over_omega_n,subs(mag,ZETA,zeta(i)),[0 omega_over_omega_n_max],...
        'DisplayName',['\zeta=' num2str(zeta(i))],'Color',(i-1)/length(zeta)*[1 1 1])

    nexttile(2)
    fplot(omega_over_omega_n,subs(phase/pi,ZETA,zeta(i)),[0 omega_over_omega_n_max],...
        'DisplayName',['\zeta=' num2str(zeta(i))],'Color',(i-1)/length(zeta)*[1 1 1])
end

%% symbolic: second order system with extra high freq zeros+poles
OMEGA_N_zeros = OMEGA_N * 1.15;
zeta_zeros = ZETA/3;
omega_n_poles = OMEGA_N * 1.2;
zeta_poles = ZETA*3;
two_zeros = (s/OMEGA_N_zeros)^2 + 2*zeta_zeros*s/OMEGA_N_zeros + 1;
two_poles = 1/( (s/omega_n_poles)^2 + 2*zeta_poles*s/omega_n_poles + 1);
X_over_F_norm_weird = X_over_F_norm * two_zeros^1.5 * two_poles^1.5;
pretty(X_over_F_norm_weird)

mag_weird = abs(X_over_F_norm_weird);
phase_weird = angle(X_over_F_norm_weird);

omega_over_omega_n_max = max(w/wn,[],'all');
for i=1:length(zeta)
    nexttile(1)
    %fplot(omega_over_omega_n,subs(mag_weird,ZETA,zeta(i)),[0 omega_over_omega_n_max],...
    %    'DisplayName',['\zeta=' num2str(zeta(i))],'Color',(i-1)/length(zeta)*[1 0 0])

    nexttile(2)
    %fplot(omega_over_omega_n,subs(phase_weird/pi,ZETA,zeta(i)),[0 omega_over_omega_n_max],...
    %    'DisplayName',['\zeta=' num2str(zeta(i))],'Color',(i-1)/length(zeta)*[1 0 0])
end

% axis labels
nexttile(1)
ylabel(['Magnitude |' y_lab '| ' y_lab_extra ' (-)'])
xlim([.5 omega_over_omega_n_max])
nexttile(2)
xlabel('\omega/\omega_n (-)')
ylabel(['Phase \angle(' y_lab ') / \pi'])
xlim([.5 omega_over_omega_n_max])
ylim([-1 0])
plot(NaN,NaN,'k*-','DisplayName','RM3') % dummy for legend

% legend, colorbar, and figure size
improvePlot
leg = legend();
leg.Layout.Tile = 'east';
f.Position(2:4) = [42 725 836];
colormap( [linspace(0,1,length(p.Hs)).', zeros(length(p.Hs),1), linspace(1,0,length(p.Hs)).'] )
cb = colorbar;
cb.Layout.Tile = 'east';
clim([p.Hs(1)-.25; p.Hs(end)+.25])
cb.Ticks = p.Hs;
cb.TickLabels = num2str(p.Hs);%.25+[p.Hs(1)-.5; p.Hs]);
cb.Label.String = 'H_s (m)';

%% estimation
mag_data   = mag_matrix(5,:);
angle_data = angle_matrix(5,:);
TF_data = mag_data .* exp(1i * angle_data);
TF_w = w(5,:);

% fill in holes with other Hs
mag_nan   = mean(mag_matrix(:,isnan(TF_data)),'omitnan');
angle_nan = mean(angle_matrix(:,isnan(TF_data)),'omitnan');
TF_data(isnan(TF_data)) = mag_nan .* exp(1i * angle_nan);
unique_w = unique(w(:,isnan(TF_w)),'stable');
TF_w(isnan(TF_w)) = unique_w(~isnan(unique_w));

TF_frd = frd(TF_data(~isnan(TF_data)), TF_w(~isnan(TF_w)));
%np = 4; nz = 2; % okay fit, stable
%np = 5; nz = 3; % good fit, unstable
for np=4:6
    np
    est_TF = tfest(TF_frd,np)%,nz);
    [Z,P,K] = zpkdata(est_TF);
    zeros = Z{:}
    poles = P{:}
    gain  = K
    stable = allmargin(est_TF).Stable
    
    figure(2)
    bode(est_TF)
    hold on
    
    figure(3)
    pzplot(est_TF)
    hold on
end
figure(2)
bode(TF_frd,'k*')
% crtlpref > response > phase wrap at -360 (could also do with
% PhaseWrappingEnabled and PhaseWrappingBranch of bodeplot command)
legend(num2str((4:6).'))
improvePlot

