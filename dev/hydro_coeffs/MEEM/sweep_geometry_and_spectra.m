clear all
close all

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

    
    % force
    force_2_over_rho_w_eiwt = 1i * int( int(r*(phi_h_i2+phi_p_i2),a1,a2), 0, 2*pi );
    
    %symvar(force_2_over_rho_w_eiwt)

    force = real(double(vpa(subs(force_2_over_rho_w_eiwt,{h,a1,a2,d2},{h_num,a1_num,a2_num,d2_num}))));
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