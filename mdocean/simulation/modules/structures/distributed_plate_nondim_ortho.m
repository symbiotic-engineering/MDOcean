function [w_nondim,Mr_nondim,Mt_nondim,Q_nondim] = distributed_plate_nondim_ortho(a,b,q,nu,rho, ...
                                                                            Dtheta_over_Dr, Drt_over_Dr)
% Roark's table 11.2, case 2L (page 467)
% a: outer radius
% b: inner radius
% r0: radius at which the distributed loading begins
% Now extended for orthotropic bending rigidity by allowing
% different radial/tangential equivalent thicknesses.
%

    % if orthotropic constants aren't given, assume isotropic (all axes the same)
    if nargin < 6
        Dtheta_over_Dr = 1;
    end
    % approximation for coupling: D_rt = nu * sqrt(Dr * Dtheta), reduces to nu*D for isotropic
    if nargin < 7
        Drt_over_Dr = nu*sqrt(Dtheta_over_Dr);
    end

    v = nu;
    r = rho*a;
    r0 = b;

    C8 = C_function(8, a, b, v);
    C9 = C_function(9, a, b, v);

    F2 = F_function(2, r, b, v);
    F3 = F_function(3, r, b, v);
    F5 = F_function(5, r, b, v);
    F6 = F_function(6, r, b, v);
    F8 = F_function(8, r, b, v);
    F9 = F_function(9, r, b, v);

    L17 = L_function(17, r0, a, Dtheta_over_Dr);

    G11 = G_function(11, r0, r, Dtheta_over_Dr);
    G14 = G_function(14, r0, r, Dtheta_over_Dr);
    G17 = G_function(17, r0, r, Dtheta_over_Dr);

    % Mrb: reaction moment per unit length on inner edge
    Mrb = -q*a^2/C8 * (C9*(a^2-r0^2)/(2*a*b) - L17);
    % Qb: shear reaction force per unit length on inner edge
    Qb  =  q/2/b * (a^2 - r0^2);

    % =================================================================
    % === ORTHOTROPIC MODIFY: nondimensional rigidity ratios ==========
    % =================================================================
    Deff_over_Dr   = Dtheta_over_Dr - Drt_over_Dr.^2;   % effective rigidity ratio

    % =================================================================
    % === y*D  and theta*D (only q-terms need D_eff) ==================
    % =================================================================

    % case 2 general deflection equation p463, with first two terms
    % ommitted since yb and theta_b are zero for case 2L
    y_times_Dr = Mrb * r.^2 .* F2 ...
               + Qb  * r.^3 .* F3 ...
               - q   * r.^4 .* G11;

    % theta_times_D:
    %   Mrb and Qb contributions use Dr
    %   q term uses D_eff
    theta_times_Dr = Mrb * r    .* F5 ...
                   + Qb  * r.^2 .* F6 ...
                   - q   * r.^3 .* (Deff_over_Dr * G14);

    % =================================================================
    % === Mr and Mt ===================================================
    % =================================================================
    % Mr depends only on radial rigidity Dr, so unchanged:
    Mr = Mrb * F8 ...
       + Qb * r .* F9 ...
       - q  * r.^2 .* G17;

    % Orthotropic tangential moment:
    % Mt = (Drt/Dr)*Mr  - (D_eff)*(theta/r)
    Mt = Drt_over_Dr .* Mr  - Deff_over_Dr .* (theta_times_Dr ./ r);

    Q = Qb * b./r ...
        - q./(2*r).*(r.^2 - r0.^2) .* bracket;

    % Roark's p463 nondimensionalization
    K_y = y_times_D / (q * a^4);
    K_Mt = Mt / (q * a^2);
    K_Mr = Mr / (q * a^2);
    K_theta = theta_times_D / (q * a^3);
    K_Q = Q / (q * a);

    % nondimensionalize to be consistent with concentrated loading plate in Boedo and Prantil
    % note: w_nondim = 2*Ky/(1-(b/a)^2)
    P = q * pi * (a^2-b^2); % equivalent concentrated load
    w_nondim  = y_times_Dr * 2*pi/(P*a^2);
    Mr_nondim = Mr * 2*pi/P;
    Mt_nondim = Mt * 2*pi/P;
    Q_nondim  = Q * 2*pi/(P*a);

    % use sign convention that positive F_heave results in positive
    % deflection at outer end and positive moment at inner end
    w_nondim  = -w_nondim;
    Mr_nondim = -Mr_nondim;
    Mt_nondim = -Mt_nondim;
    Q_nondim  = -Q_nondim;
end

function L = L_function(number, r0, a, Dtheta_over_Dr)
    ratio = r0/a;
    L = LG_function(number, ratio, Dtheta_over_Dr);
end

function G = G_function(number, r0, r, Dtheta_over_Dr)

    bracket = zeros(size(r));
    bracket(r > r0) = 1;

    ratio = r0./r;

    G = LG_function(number, ratio, Dtheta_over_Dr) .* bracket;

end

function F = F_function(number, r, b, v)
    ratio = r/b;
    F = CF_function(number, ratio, v);
end

function C = C_function(number, a, b, v)
    ratio = a/b;
    C = CF_function(number, ratio, v);
end

function CF = CF_function(number, ratio, v)

    if number==2
        CF = 1/4 * (1 - (1./ratio).^2.*(1+2*log(ratio)));
    elseif number==3
        CF = 1./ratio/(4) .* ( ( (1./ratio).^2 + 1 ).*log(ratio) + (1./ratio).^2 - 1);
    elseif number==5
        CF = 1/2 * (1 - (1./ratio).^2);
    elseif number==6
        CF = 1./ratio/(4) .* ( (1./ratio).^2 - 1 + 2*log(ratio) );
    elseif number==8
        CF = 1/2 * (1+v+(1-v)*(1./ratio).^2);
    elseif number==9
        CF = 1./ratio .* (1/2*(1+v)*(log(ratio)) + 1/4*(1-v)*(1-(1./ratio).^2));
    else
        error(['number=' number ' is not implemented'])
    end
end
function LG = LG_function(number, ratio, Dtheta_over_Dr)

    if number==11
        LG = 1/(8*(9-Dtheta_over_Dr)) * (1 + 4*ratio.^2 - 5*ratio.^4 - 4*ratio.^2.*(2+ratio.^2).*log(1./ratio) );
    elseif number==14
        LG = 1/(2*(9-Dtheta_over_Dr)) * (1 - (ratio).^4  - 4*(ratio).^2 .* log(1./ratio));
    elseif number==17
        LG = 2/(9-Dtheta_over_Dr) * (1 - ((1-v)/4)*(1-ratio.^4) - (ratio).^2.*(1+(1+v)*log(1./ratio)));
    else
        error(['number=' number ' is not implemented'])
    end

end
