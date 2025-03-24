clear
close all

% sweep variables
R = 0.1 : 0.02 : 10; % normalized resistance ratio
D = 0 : 0.2 : 1; % ratio of hydrodynamic damping to total damping
alpha_m = -10 : 0.1 : 10; % reactive to real impedance ratio
% F_max_over_Fp = 0.3 : 0.2 : 1.1; % sets level of saturation

% generate sweep grid
[R_grid, D_grid, ALPHA_M_grid] = meshgrid(R, D, alpha_m);

% constant variables - changing these should have no effect on the nondimensional outputs
m = 1;  % dimensional term for mass
w = 1;  % dimensional term for frequency
F_h = 1; % dimensional term for exciting force

% get results
tic
[avg_pwr, voltage_ratio, current_ratio] = describing_fcn(R_grid, D_grid,ALPHA_M_grid, 7);
toc
%tic
% [avg_pwr2, max_x2, max_xdot2, pwr_ratio2, x_ratio2, xdot_ratio2] = ground_truth_wot(  ZETA_U, W_U_STAR, F_MAX_FP, m, w, F_h);
%[avg_pwr2, max_x2, max_xdot2, pwr_ratio2, x_ratio2, xdot_ratio2] = ground_truth(  ZETA_U, W_U_STAR, F_MAX_FP, m, w, F_h);
%toc

% combine results
dim_cat = ndims(R_grid)+1;
results_sim = cat(dim_cat, avg_pwr, voltage_ratio, current_ratio);
% results_act = cat(dim_cat, avg_pwr2, max_x2, max_xdot2, pwr_ratio2, x_ratio2, xdot_ratio2);
results_str = {'Average Power', 'Voltage Ratio', 'Currrent Ratio'};

%results_unsat_sim = results_sim(F_MAX_FP >= 1);
% results_unsat_act = results_act(F_MAX_FP >= 1);

for j = 1:size(results_sim,dim_cat)
    figure
    res_sim = results_sim(:,:,:,j);
    % res_act = results_act(:,:,:,j);
    % color_lims = [min([res_sim,res_act],[],'all') max([res_sim,res_act],[],'all')];

    % plot simulated results
    subplot(1,1,1)
    slice(R_grid, D_grid, ALPHA_M_grid, res_sim,...
        R(2:end),D(end),alpha_m);
    % caxis(color_lims)
    shading interp
    colorbar

    xlabel('R');
    ylabel('D');
    zlabel('\alpha_m');
    title('Describing Function')
    axis tight
    grid on
    
    % plot actual results
    % subplot 132
    % slice(ZETA_U,W_U_STAR,F_MAX_FP, res_act,...
    %     zeta_u(2:end),w_u_star(end),F_max_over_Fp,'nearest');
    % caxis(color_lims)
    % colorbar
    % xlabel('\zeta_u')
    % ylabel('\omega_u^*')
    % zlabel('F_{max}/F_p')
    % title('Numerical Integration')
    % 
    % error_h = subplot(1,3,3);
    % slice(ZETA_U,W_U_STAR,F_MAX_FP, (res_sim-res_act)./res_act,...
    %     zeta_u(2:end),w_u_star(end),F_max_over_Fp,'nearest');
    % caxis([-1,1])
    % colormap(error_h,bluewhitered)
    % colorbar
    % xlabel('\zeta_u')
    % ylabel('\omega_u^*')
    % zlabel('F_{max}/F_p')
    % title('Fractional Error')

    sgtitle(results_str{j})
end

function [f_sat_i] = get_f_sat(F_max_over_Fp, i)
    assert(all(F_max_over_Fp <= 1 & F_max_over_Fp >= 0,'all'))
    if i == 1
        f_sat_i = 2/pi * (F_max_over_Fp .* sqrt(1 - F_max_over_Fp.^2) + asin(F_max_over_Fp));

    elseif rem(i,2)==1 && i>0 % odd positive numbers
        theta = i * asin(F_max_over_Fp);
        f_sat_i = 4/pi * 1/(i*(i^2-1)) * (i * sqrt(1 - F_max_over_Fp.^2) .* sin(theta) ...
            - F_max_over_Fp .* cos(theta));

    else % nonpositive numbers and even positive numbers
        f_sat_i = 0;
    end

end

function [avg_pwr, voltage_ratio, current_ratio] = describing_fcn(R, D, ...
                                                alpha_m, order)

    % Compute unsaturated power and voltage ratio based on Force_saturation_conference equations
    G_0 = 1; % constant from the paper
    J = 1;  % placeholder for incident energy density
    k = 1;  % placeholder for wavenumber
    L=0;
    F_e=0;
    G=1;
    Kt=1;

    Z_th = R .* ( 1 + 1i*L ) + D ./ ( R .* (1 + 1i*alpha_m) );
    V_th = sqrt(F_e.*R./(Kt.*G)).*(D./(R.*sqrt(1 + alpha_m.^2)));

    z=Z_th;

    Gamma = (z+1)./(z-1); % Define or compute Gamma as needed
    

    % Calculate Re(Gamma) and Im(Gamma)
    Re_Gamma = real(Gamma);
    Im_Gamma = imag(Gamma);

    alpha = imag(z)./real(z);

    epsilon = 1; % Define epsilon
    % Calculate voltage_ratio based on the given formula
    voltage_ratio = sqrt((abs(Gamma).^2 + 2 .* epsilon .* Re_Gamma + 1) ./ ...
                           (alpha.^2 .* abs(Gamma).^2 + 2 .* alpha .* Im_Gamma + 1));

    epsilon = -1; % Define epsilon
    % Calculate current_ratio based on the given formula
    current_ratio = sqrt((abs(Gamma).^2 + 2 .* epsilon .* Re_Gamma + 1) ./ ...
                           (alpha.^2 .* abs(Gamma).^2 + 2 * alpha .* Im_Gamma + 1));

    % P_unsat corresponds to P_L^m in matched case (Equation 7)
    P_unsat = (G_0 .* J ./ k) * (D ./ (1 + (R ./ D) .* (1 + alpha_m.^2)));
    P_L = P_unsat .* (1-abs(Gamma).^2);

    avg_pwr = P_L;


end



