clear all
close all

N = 20;
syms r theta z real                                     % coordinates
syms R_1n_1(n) R_1n_2(n) R_2n_2(n) Z_n_i1(n)  ...       % basis functions
        Z_n_i2(n) Lambda_k(k) N_k(k) Z_k_e(k) m_k(k)              
syms h a1 a2 d1 d2 m0 real positive                     % constants
syms n k real positive                                  % counter indices
syms C_1n_1(n) C_1n_2(n) C_2n_1(n) C_2n_2(n) B_k(k)     % unknown coefficients

% equation numbers refer to Chau & Yeung 2010 unless otherwise noted

%% setup analytical boundary value problem equations
% eq 4
lambda_n1(n) = n*pi/(h-d1);
lambda_n2(n) = n*pi/(h-d2);

% eq 5
phi_p_i1 = 1/(2*(h-d1)) * ((z+h)^2 - r^2/2);
phi_p_i2 = 1/(2*(h-d2)) * ((z+h)^2 - r^2/2);

% eq 7
R_1n_1(n) = piecewise(n==0, 1/2, n>=1, ...
                    besseli(0,lambda_n1(n)*r)/besseli(0,lambda_n1(n)*a2));
R_1n_2(n) = piecewise(n==0, 1/2, n>=1, ...
                    besseli(0,lambda_n2(n)*r)/besseli(0,lambda_n2(n)*a2));

% eq 8
R_2n_1(n) = sym(0);
R_2n_2(n) = piecewise(n==0, 1/2*log(r/a2), ...
    n>=1, besselk(0,lambda_n2(n)*r)/besselk(0,lambda_n2(n)*a2));

% eq 9
Z_n_i1(n) = piecewise(n==0, 1, n>=1, sqrt(2)*cos(lambda_n1(n)*(z+h)));
Z_n_i2(n) = piecewise(n==0, 1, n>=1, sqrt(2)*cos(lambda_n2(n)*(z+h)));

% eq 6
phi_h_n_i1(n) = (C_1n_1(n) * R_1n_1(n) + C_2n_1(n) * R_2n_1(n)) * Z_n_i1(n);
phi_h_n_i2(n) = (C_1n_2(n) * R_1n_2(n) + C_2n_2(n) * R_2n_2(n)) * Z_n_i2(n);

% eq 13
Lambda_k(k) = piecewise(k==0, besselh(0,1,m0*r)/besselh(0,1,m0*a2), ...
    k>=1, besselk(0,m_k(k)*r)/besselk(0,m_k(k)*a2));

% eq 2.34 in analytical methods book, also eq 16 in Seah and Yeung 2006 
N_k(k) = piecewise(k==0, 1/2*(1+sinh(2*m0*h)/(2*m0*h)), ...
                   k>=1, 1/2*(1+sin(2*m_k(k)*h)/(2*m_k(k)*h)) );

% eq 14
Z_k_e(k) = piecewise(k==0, 1/sqrt(N_k(k)) * cosh(m0 * (z+h)), ...
                     k>=1, 1/sqrt(N_k(k)) * cos(m_k(k) * (z+h)));

B_n(n) = subs(B_k, k, n);
Lambda_n(n) = subs(Lambda_k, k, n);
Z_n_e(n) = subs(Z_k_e, k, n);

% potential matching
% equation 22 in old 1981 paper, applied to boundary 2-e
dz_2 = h - d2;
match_2e_potential = (C_1n_2(n) * subs(R_1n_2(n),r,a2) + ...
                     C_2n_2(n) * subs(R_2n_2(n),r,a2)) / dz_2 == ...
    B_n(n) * subs(Lambda_n(n),r,a2) / dz_2 - ...
    int(subs(phi_p_i2,r,a2) * Z_n_i2(n), z, -h, -d2);

% equation 22 in old 1981 paper, applied to boundary 1-2
dz_1 = h - d1;
match_12_potential = C_1n_1(n) * subs(R_1n_1(n),r,a1) / dz_1 == ...
    ( C_1n_2(n) * subs(R_1n_2(n),r,a1) + ...
      C_2n_2(n) * subs(R_2n_2(n),r,a1) ) / dz_1 + ...
    int(subs(phi_p_i2 - phi_p_i1,r,a1) * Z_n_i1(n), z, -h, -d1);

% velocity matching
% equation 23 in old 1981 paper, applied to boundary 2-e
match_2e_velocity = B_n(n) * subs(diff(Lambda_n(n), r), r, a2) / dz_2 == ...
    (C_1n_2(n) * subs(diff(R_1n_2(n), r), r, a2) + ...
     C_2n_2(n) * subs(diff(R_2n_2(n), r), r, a2) ) / dz_2 + ...
    int(subs(diff(phi_p_i2,r),r,a2) * Z_n_e(n), z, -h, -d2 );

% equation 23 in old 1981 paper, applied to boundary 1-2
match_12_velocity = (C_1n_2(n) * subs(diff(R_1n_2(n),r), r, a1) + ...
    C_2n_2(n) * subs(diff(R_2n_2(n),r), r, a1)) / dz_1 == ...
    C_1n_1(n) * subs( diff(R_1n_1(n),r), r,a1) / dz_1 + ...
    int( subs(diff(phi_p_i1 - phi_p_i2,r),r,a1) * Z_n_i2(n), z, -h, -d1 );

% substitute in 0:N for count variables
eqns = [subs(match_12_potential,n,0:N), subs(match_2e_potential,n,0:N) ...
        subs(match_12_velocity,n,0:N),  subs(match_2e_velocity,n,0:N)];

Lambda_k = subs(Lambda_k,k,0:N);
Z_k_e = subs(Z_k_e,k,0:N);
unknowns = [C_1n_1(0:N) C_1n_2(0:N) C_2n_2(0:N) B_n(0:N)];

% substitute symbolic constant vectors for symbolic functions
syms C_1n_1_const C_1n_2_const C_2n_2_const B_k_const [N+1 1] real
unknowns_const = [C_1n_1_const; C_1n_2_const; C_2n_2_const; B_k_const];
eqns = subs(eqns, unknowns, unknowns_const');

syms m_k_const [1 N]
eqns = subs(eqns,m_k(1:N),m_k_const);
Lambda_k = subs(Lambda_k,m_k(1:N),m_k_const);
Z_k_e = subs(Z_k_e,m_k(1:N),m_k_const);

%symvar(eqns)

%% prepare geometries to numerical sweep
h_mat = 10;

spatial_res = 30;
sweep_res = 10;

a2_vec = h_mat*1; % linspace(.5, 1, sweep_res);
d2_vec = h_mat*[.25];% .4];
a1_vec = h_mat*.5;
d1_vec = h_mat*[.5];% .75];
m0_vec = [.1];% 2];

% check valid geometry
assert(max(d1_vec) < min(h_mat) && ...
       max(d2_vec) < min(h_mat) && ...
       max(a1_vec) < min(a2_vec) && ...
       max(d2_vec) < min(d1_vec) )

[a2_mat,d2_mat,a1_mat,d1_mat] = ndgrid(a2_vec,d2_vec,a1_vec, d1_vec);

phi = zeros(spatial_res,spatial_res+2,numel(a2_mat),length(m0_vec));
[force, A, B, sigma_force, sigma_position] = deal(zeros(numel(a2_mat), length(m0_vec)));

plot_phi = true;
plot_omega = true;

% eq 11 - numerical root-finding to solve the imaginary dispersion relation
%figure
%hold on
m_k_mat = zeros(length(m0_vec),N);

for freq_idx = 1:length(m0_vec)
    m_k_vec = zeros(1,N);
    m0_i = m0_vec(freq_idx);
    m_k_sq_err = @(m_k) (m_k * tan(m_k * h_mat) + m0_i * tanh(m0_i*h_mat))^2;
    for k_idx=1:N
        m_k_lower = (pi*(k_idx-1) + pi/2) / m0_i + eps;
        m_k_upper = (pi*(k_idx-1) + pi)   / m0_i - eps;
        m_k_vec(k_idx) = fminbnd(m_k_sq_err,m_k_lower,m_k_upper);
    end
    %plot(1:N,m_k_vec,'*-')
    %xlabel('k'); ylabel('m_k')
    shouldnt_be_int = round(m0_i * m_k_vec/pi - .5, 4); % check in case it found an asymptote instead of a solution
    not_repeated = numel(unique(m_k_vec)) == numel(m_k_vec);
    assert(all(shouldnt_be_int ~= floor(shouldnt_be_int)) && not_repeated)

    m_k_mat(freq_idx,:) = m_k_vec;
end

%% numerically solve equations for geometry sweep
for geom_idx = 1:numel(a2_mat)
    disp(['running geometry ' num2str(geom_idx) ' out of ' num2str(numel(a2_mat))]);
    
    for freq_idx = 1:length(m0_vec)
        disp(['    running frequency ' num2str(freq_idx) ' out of ' num2str(length(m0_vec))]);

        [phi(:,:,geom_idx,freq_idx), force(geom_idx,freq_idx) ] = get_phi_force(...
                                                    eqns, unknowns_const, N, ...
            R_1n_1, R_1n_2, R_2n_1, R_2n_2, Z_n_i1, Z_n_i2, Lambda_k, Z_k_e, ...
            phi_p_i1, phi_p_i2, r, z, spatial_res, plot_phi, ...
            h_mat, m_k_mat(freq_idx,:), a1_mat(geom_idx), a2_mat(geom_idx), ...
            d1_mat(geom_idx), d2_mat(geom_idx), m0_vec(freq_idx));
    end

    [A(geom_idx,:), B(geom_idx,:), sigma_force(geom_idx,:), sigma_position(geom_idx,:)] = ...
        get_frequency_stuff(force(geom_idx,:), m0_vec, h_mat, plot_omega);

end

%% multidim plot over geometry
if length(d2_vec)>1 && length(d1_vec)>1
    figure
    contourf(reshape(d1_mat,length(d2_vec),length(d1_vec)),...
             reshape(d2_mat,length(d2_vec),length(d1_vec)),...
             reshape(force(:,1),length(d2_vec),length(d1_vec)))
    xlabel('d_1')
    ylabel('d_2')
    colorbar
    title('Heave Exciting Force: X/\rho \omega e^{i\omega t}')
end

%%
function [phi, force] = get_phi_force(eqns, unknowns_const, N, ...
                R_1n_1, R_1n_2, R_2n_1, R_2n_2, Z_n_i1, Z_n_i2, Lambda_k, Z_k_e, ...
                phi_p_i1, phi_p_i2, r, z, spatial_res, plot_phi, ...
                h_num, m_k_num, a1_num, a2_num, d1_num, d2_num, m0_num)

    syms m_k_const [1 N]
    eqns = subs(eqns,m_k_const,m_k_num);

    syms h a1 a2 d1 d2 m0 real positive 
    params = {h a1 a2 d1 d2 m0};
    params_num = {h_num a1_num a2_num d1_num d2_num m0_num};
    eqns = subs(eqns,params,params_num);
    
    solns = vpasolve(eqns, unknowns_const);
    
    fns = fieldnames(solns);
    
    C_1n_1s = ones(N+1,1);
    C_2n_1s = zeros(N+1,1); % because R_2n_1s = 0
    C_1n_2s = ones(N+1,1);
    C_2n_2s = ones(N+1,1);
    B_ks = ones(N+1,1);
    index = 1;
    
    for k = 1:length(fns)
        if index > N+1
            index = 1;
        end
        if isempty(solns.(fns{k}))
            error('Linear solver could not find a solution')
        else
            if (k<=N+1)
                C_1n_1s(index) = solns.(fns{k});
            end
            if (N+1<k) && (k<=2*(N + 1))
                C_1n_2s(index) = solns.(fns{k});   
            end
            if (k>2*(N+1)) && (k<=3*(N+1))
                 C_2n_2s(index) = solns.(fns{k});      
            end
            if (k>3*(N+1)) 
                 B_ks(index) = solns.(fns{k});
            end
        end
        index = index + 1;   
    end
    
    phi_h_n_i1_solns_all = (C_1n_1s' .* R_1n_1(0:N) + C_2n_1s' .* R_2n_1(0:N)) .* Z_n_i1(0:N); 
    phi_h_n_i2_solns_all = (C_1n_2s' .* R_1n_2(0:N) + C_2n_2s' .* R_2n_2(0:N)) .* Z_n_i2(0:N);
    
    phi_e_k = B_ks' .* Lambda_k .* Z_k_e;
    
    % summing all to get the potential 
    phi_h_i1 = sum(phi_h_n_i1_solns_all,2);
    phi_h_i2 = sum(phi_h_n_i2_solns_all,2);
    
    phi_e_sym = sum(phi_e_k);
    
    % sub in geometry
    phi_h1 = subs(phi_h_i1,params,params_num);
    phi_h2 = subs(phi_h_i2,params,params_num);
    phi_p1 = subs(phi_p_i1,params,params_num);
    phi_p2 = subs(phi_p_i2,params,params_num);
    phi_e = subs(phi_e_sym,params,params_num);
    phi_e = subs(phi_e,m_k_const,m_k_num);

    % 36 and 37
    phi_1 = phi_h1 + phi_p1;
    phi_2 = phi_h2 + phi_p2;
    
    % velocity
    v_1r = diff(phi_1,r);
    v_1z = diff(phi_1,z);
    v_2r = diff(phi_2,r);
    v_2z = diff(phi_2,z);
    v_er = diff(phi_e,r);
    v_ez = diff(phi_e,z);

    % sub in spatial coordinates
    r_vec = linspace(2*a2_num/spatial_res,2*a2_num,spatial_res);
    r_vec = sort([r_vec a1_num a2_num]);
    z_vec = linspace(-h_num,0,spatial_res);
    [R,Z] = meshgrid(r_vec,z_vec);
    
    phi1 = double(subs(phi_1,{r,z},{R,Z}));
    phi2 = double(subs(phi_2,{r,z},{R,Z}));
    phi1h = double(subs(phi_h1,{r,z},{R,Z}));
    phi1p = double(subs(phi_p1,{r,z},{R,Z}));
    phi2h = double(subs(phi_h2,{r,z},{R,Z}));
    phi2p = double(subs(phi_p2,{r,z},{R,Z}));
    phie = double(subs(phi_e,{r,z},{R,Z}));

    v1r = double(subs(v_1r,{r,z},{R,Z}));
    v1z = double(subs(v_1z,{r,z},{R,Z}));
    v2r = double(subs(v_2r,{r,z},{R,Z}));
    v2z = double(subs(v_2z,{r,z},{R,Z}));
    ver = double(subs(v_er,{r,z},{R,Z}));
    vez = double(subs(v_ez,{r,z},{R,Z}));
    
    % assemble total phi based on phi in each region
    regione = R > a2_num;
    region1 = R <= a1_num & Z < -d1_num;
    region2 = R > a1_num & R <= a2_num & Z < -d2_num;
    
    phi = NaN(size(R));
    phi(region1) = phi1(region1);
    phi(region2) = phi2(region2);
    phi(regione) = phie(regione);

    phiH = NaN(size(R));
    phiH(region1) = phi1h(region1);
    phiH(region2) = phi2h(region2);

    phiP = NaN(size(R));
    phiP(region1) = phi1p(region1);
    phiP(region2) = phi2p(region2);

    v_r = NaN(size(R));
    v_r(region1) = v1r(region1);
    v_r(region2) = v2r(region2);
    v_r(regione) = ver(regione);

    v_z = NaN(size(R));
    v_z(region1) = v1z(region1);
    v_z(region2) = v2z(region2);
    v_z(regione) = vez(regione);

    if plot_phi
        region_body = ~region1 & ~region2 & ~regione;
        plot_potential(phi,R,Z,region_body,'Total Potential');
        plot_potential(phiH,R,Z,region_body,'Homogeneous Potential');
        plot_potential(phiP,R,Z,region_body,'Particular Potential');
        plot_potential(v_r,R,Z,region_body,'Radial Velocity')
        plot_potential(v_z,R,Z,region_body,'Vertical Velocity')
        plot_velocity(real(v_r),real(v_z),R,Z);

        plot_matching(phi1,phi2,phie,a1_num,a2_num,R,Z,'\phi')
        plot_matching(v1r,v2r,ver,a1_num,a2_num,R,Z,'v_r')
    end

    % force
    force_2_over_rho_w_eiwt = 1i * int( int(r*(phi_h_i2+phi_p_i2),a1,a2), 0, 2*pi );
    
    %symvar(force_2_over_rho_w_eiwt)

    force = real(double(vpa(subs(force_2_over_rho_w_eiwt,{h,a1,a2,d2},{h_num,a1_num,a2_num,d2_num}))));
end

function [] = plot_potential(phi,R,Z,region_body,name)
    figure
    minphi = min(real(phi),[],'all');
    maxphi = max(real(phi),[],'all');

    num_levels = 30;
    if minphi<0
        levels = [-logspace(log10(-minphi),0,num_levels/2), logspace(0,log10(maxphi),num_levels/2)];
        levels = sort([levels, 100:5:180]);
    else
        levels = linspace(minphi,maxphi,num_levels);%[logspace(log10(minphi),log10(maxphi),num_levels)];
    end
    subplot 121
    [c,h_fig] = contourf(R,Z,real(phi),levels);
    clabel(c,h_fig)
    xlabel('R')
    ylabel('Z')
    title([name ' - Real'])
    colorbar;
    set(gca,'ColorScale','log')
    
    imag_phi = imag(phi);
    if numel(unique(imag_phi)) > 1
        imag_phi(region_body) = NaN;
        subplot 122
        [c,h_fig] = contourf(R,Z,imag_phi);
        clabel(c,h_fig)
        xlabel('R')
        ylabel('Z')
        title([name ' - Imaginary'])
        colorbar
    end
end
function plot_velocity(v_r,v_z,R,Z)
    v_tot = sqrt(v_r.^2 + v_z.^2);

    num_levels = 10;
    levels = logspace(1,6,num_levels);%logspace(log10(min(v_tot,[],'all')),log10(max(v_tot,[],'all')),num_levels);
    levels = sort([levels, 0:10, 12:2:20]);

    figure
    [c,h_fig] = contourf(R,Z,v_tot,levels);
    clabel(c,h_fig)
    xlabel('R')
    ylabel('Z')
    title('Velocity')
    colorbar
    hold on
    quiver(R,Z,v_r./v_tot,v_z./v_tot)
    set(gca,'ColorScale','log')
end

function plot_matching(phi1,phi2,phie,a1,a2,R,Z,name)
    % at R = a1
    [~,idx_a1] = min(abs(R - a1),[],2,'linear');
    phi1_a1 = abs(phi1(idx_a1));
    phi2_a1 = abs(phi2(idx_a1));

    % at R = a2
    [~,idx_a2] = min(abs(R - a2),[],2,'linear');
    phi2_a2 = abs(phi2(idx_a2));
    phie_a2 = abs(phie(idx_a2));

    % plot
    figure
    plot(Z(idx_a1),phi1_a1,'r--',Z(idx_a1),phi2_a1,'m-')
    hold on
    plot(Z(idx_a2),phi2_a2,'b-',Z(idx_a2),phie_a2,'c--')
    legend('phi_1,a1','phi_2,a1','phi_2,a2','phi_e,a2')
    xlabel('Z')
    ylabel(['|' name '|'])
    title([name ' Matching'])
end

function [A,B,sigma_force,sigma_position] = get_frequency_stuff(force,k,h,plot_omega)
    g = 9.8;
    rho = 1025;
    omega = sqrt(k*g.*tanh(k*h));
    X3 = force * rho .* omega;
    Vg = g./(2*omega) .* (k.^2 + tanh(k*h) - k.^2.*tanh(k*h).^2);
    B = k./ (4*rho*g*Vg) .* X3.^2;

    A = zeros(size(B));
    for i=1:length(omega)
        A(i) = 1/omega(i)^2 * hilbert(-omega(i)*B(i)); % kramers kronig relationship
    end

    if length(omega) > 1
        Hm0 = 2;
        Tp = 1;
        amplitude = Hm0 / 2;
        H_force = X3 / amplitude;
    
        % transfer function
        m = 1; % mass
        A = 1;
        C = rho * g * A;
        H_position = X3 ./ (-omega.^2.*(m+A) + C + 1i*omega.*B);
        if plot_omega
            figure
            plot(omega,abs(H_position))
            title('Heave RAO')
            xlabel('\omega (rad/s)')
            ylabel('|\Xi/A| - Magnitude of Position Transfer Function')
        end
    
        % spectrum and standard deviations
        S = pierson(omega, Hm0, Tp);                    % spectrum
        sigma_force = wkinchin(S,omega,H_force);        % standard deviation force
        sigma_position = wkinchin(S,omega,H_position);  % standard deviation pos
    else
        sigma_force = 0;
        sigma_position = 0;
    end
end

function E_pm = pierson(omega,Hm0, Tp)
    gamma = 3.3;
    alpha  = 1. / (.23 + .03 * gamma - .185 / (1.9 + gamma)) / 16. ; 
    E_pm = alpha .* Hm0^2 * Tp^-4 .* omega.^-5 .* exp(-1.25 * (Tp .* omega).^-4);
end

function std = wkinchin(S, omega, H)
    var = trapz( omega, S .* abs(H).^2 );
    std = sqrt(var);
end