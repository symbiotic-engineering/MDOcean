clear
close all

% sweep variables
R = 0.1 : 0.02 : 10; % normalized resistance ratio
D = 0 : 0.2 : 1; % ratio of hydrodynamic damping to total damping
alpha_m = -10 : 0.1 : 10; % reactive to real impedance ratio
F_max_over_Fp = 0.3 : 0.2 : 1; % sets level of saturation

% generate sweep grid
[R_grid, D_grid, ALPHA_M_grid,F_MAX_FP] = ndgrid(R, D, alpha_m, F_max_over_Fp);

[R_grid, D_grid, ALPHA_M_grid,F_MAX_FP] = ndgrid_to_meshgrid(R_grid, D_grid, ALPHA_M_grid,F_MAX_FP);

% constant variables - changing these should have no effect on the nondimensional outputs
m = 1;  % dimensional term for mass
w = 1;  % dimensional term for frequency
F_h = 1; % dimensional term for exciting force

% Get results
tic
[avg_pwr, voltage_ratio, current_ratio] = describing_fcn(R_grid, D_grid, ALPHA_M_grid, F_MAX_FP, 7);
toc
%tic
%[avg_pwr2, max_voltage_ratio2, max_xdot2, pwr_ratio2, voltage_ratio2, xdot_ratio2] = ground_truth_wot(R_grid, D_grid, ALPHA_M_grid, F_MAX_FP, m, w, F_h);
%toc

% Combine results
dim_cat = ndims(R_grid) + 1;
results_sim = cat(dim_cat, avg_pwr, voltage_ratio, current_ratio);
%results_act = cat(dim_cat, avg_pwr2, max_voltage_ratio2, max_xdot2, pwr_ratio2, voltage_ratio2, xdot_ratio2);
results_str = {'Average Power', 'Voltage Ratio', 'Currrent Ratio'};

% Only unsaturated results
results_unsat_sim = results_sim(F_MAX_FP >= 1);
%results_unsat_act = results_act(F_MAX_FP >= 1);


% Plot each result
for j = 1:size(results_sim, dim_cat)
    figure;
    res_sim = results_sim(:,:,:,j);

    % Simulated results plot
    subplot(1, 3, 1);

    % Correctly use ndgrid variables for slicing
    slice(squeeze(R_grid(:,:,:,1)), squeeze(D_grid(:,:,:,1)), squeeze(F_MAX_FP(:,:,:,1)), res_sim, R, D, F_max_over_Fp, 'nearest')
    colorbar;
    xlabel('R');
    ylabel('D');
    zlabel('F_{max}/F_p');
    title('Describing Function');


    slice(R_grid, D_grid, F_MAX_FP, res_sim, R, D, F_max_over_Fp, 'nearest');
    %caxis(color_lims);
    colorbar;
    xlabel('R');
    ylabel('D');
    zlabel('F_{max}/F_p');
    title('Describing Function');

    % Actual results plot
    %subplot(1, 3, 2);
    %slice(R_grid, D_grid, F_MAX_FP, res_act, R(2:end), D(end), F_max_over_Fp, 'nearest');
    %caxis(color_lims);
    %colorbar;
    %xlabel('R');
    %ylabel('D');
    %zlabel('F_{max}/F_p');
    %title('Numerical Integration');

    % Error plot
    %error_h = subplot(1, 3, 3);
    %slice(R_grid, D_grid, F_MAX_FP, (res_sim - res_act) ./ res_act, R(2:end), D(end), F_max_over_Fp, 'nearest');
    %caxis([-1, 1]);
    %colormap(error_h, bluewhitered);
    %colorbar;
    %xlabel('R');
    %ylabel('D');
    %zlabel('F_{max}/F_p');
    %title('Fractional Error');

    %sgtitle(results_str{j});
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
                                                alpha_m, F_max_over_Fp, order)

    % Compute unsaturated power and voltage ratio based on Force_saturation_conference equations
    G_0 = 1; % constant from the paper
    J = 1;  % placeholder for incident energy density
    k = 1;  % placeholder for wavenumber
    L=0;
    F_e=0;
    G=1;
    Kt=1;

    Z_th = R.*(1+L.*1i+D./(R.*(1+alpha_m.*1i)));
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
    currnent_ratio = sqrt((abs(Gamma).^2 + 2 .* epsilon .* Re_Gamma + 1) ./ ...
                           (alpha.^2 .* abs(Gamma).^2 + 2 * alpha .* Im_Gamma + 1));

    % P_unsat corresponds to P_L^m in matched case (Equation 7)
    P_unsat = (G_0 .* J ./ k) * (D ./ (1 + (R ./ D) .* (1 + alpha_m.^2)));
    P_L = P_unsat .* (1-abs(Gamma).^2);

    orders = 1 : 2 : order;
    [Power_ratio, Voltage_ratio, Current_ratio, f_sat] = deal(zeros([length(orders), size(F_max_over_Fp)]));

    for n = 1:length(orders)
        i = orders(n);
        f_sat_i = get_f_sat(F_max_over_Fp, i);


        Power_ratio(n, :) = P_L(:)./P_unsat(:);
        Voltage_ratio(n, :) = voltage_ratio(:);
        

        Current_ratio(n, :) = currnent_ratio(:);
        f_sat(n, :) = f_sat_i(:);

    end

    avg_pwr = squeeze(sum(Power_ratio));
    voltage_ratio = squeeze(Voltage_ratio(1,:,:,:,:));
    disp(size(voltage_ratio));
    current_ratio = squeeze(Current_ratio(1,:,:,:,:));
    
end

function [X_mesh, Y_mesh, Z_mesh, W_mesh] = ndgrid_to_meshgrid(X, Y, Z, W)
    X_mesh = permute(X, [2, 1, 3, 4]);
    Y_mesh = permute(Y, [2, 1, 3, 4]);
    Z_mesh = permute(Z, [2, 1, 3, 4]);
    W_mesh = permute(W, [2, 1, 3, 4]);
end

% function [avg_pwr, max_x, max_xdot, ...
%     pwr_ratio, x_ratio, xdot_ratio] = ground_truth(ZETA_U, W_U_STAR, F_MAX_FP, m, w, F_h)
% 
% [avg_pwr, max_x, max_xdot, pwr_ratio, x_ratio, xdot_ratio] = deal(zeros(size(ZETA_U)));
% 
% for i = 1:numel(ZETA_U)
%     zeta_u = ZETA_U(i);
%     w_u_star = W_U_STAR(i);
%     F_max_over_Fp = F_MAX_FP(i);
% 
%     [avg_pwr(i), max_x(i), max_xdot(i)] = run_ode(zeta_u, w_u_star, F_max_over_Fp, m, w, F_h, false);
%     [avg_pwr_unsat, max_x_unsat, max_xdot_unsat] = run_ode(zeta_u, w_u_star, 1e8, m, w, F_h, false);
% 
%     pwr_ratio(i) = avg_pwr(i) / avg_pwr_unsat;
%     x_ratio(i) = max_x(i) / max_x_unsat;
%     xdot_ratio(i) = max_xdot(i) / max_xdot_unsat;
% end

% end

% function [avg_pwr, max_x, max_xdot, ...
%     pwr_ratio, x_ratio, xdot_ratio] = ground_truth_wot(ZETA_U, W_U_STAR, F_MAX_FP, m_, w_, F_h_) %#ok<STOUT> 
% 
%     % load WOT python results from mat file 
%     vars = {'avg_pwr', 'max_x', 'max_xdot', 'pwr_ratio', 'x_ratio', ...
%         'xdot_ratio','zeta_u', 'w_u_star', 'f_max_Fp', 'm', 'w', 'F_h'};
%     load('wot_sweep_results_20240102_023554_N=11.mat',vars{:})
% 
%     % confirm that WOT was run with same inputs as matlab
%     assert(all(ismembertol(ZETA_U, zeta_u),'all'))
%     assert(all(ismembertol(W_U_STAR, w_u_star),'all'))
%     assert(all(ismembertol(F_MAX_FP, f_max_Fp),'all'))
%     assert(m == m_)
%     assert(w == w_)
%     assert(F_h == F_h_)
% 
% end

