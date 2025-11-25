function [w_nondim,Mr_nondim,Mt_nondim,Q_nondim] = distributed_plate_nondim(a,b,q,nu,rho)
% Roark's table 11.2, case 2L (page 467)
% a: outer radius
% b: inner radius
% r0: radius at which the distributed loading begins

%     A = pi*a^2;
%     q = F_heave/A;
%     P = F_heave;
    v = nu;
    r = rho*a;
    r0 = b;

    C2 = 1/4 * ( 1 - (b/a)^2 * (1 + 2*log(a/b)) );
    C3 = b/4/a*(((b/a)^2+1)*log(a/b)+(b/a)^2-1);
    C8 = 1/2 * (1+v+(1-v)*(b/a)^2);
    C9 = b/a * ( (1+v)/2*log(a/b) + (1-v)/4*(1-(b/a)^2) );
    %L3 = r0/4/a*(((r0/a)^2+1)*log(a/r0)+(r0/a)^2-1); % for case 1L
    %L9 = r0/a*((1+v)/2*log(a/r0)+(1-v)/4*(1-(r0/a)^2)); % for case 1L
    L11 = 1/64 * (1+4*(r0/a)^2-5*(r0/a)^4-4*(r0/a)^2*(2+(r0/a)^2)*log(a/r0)); % for end deflection only
    L17 = 1/4 * (1-(1-v)/4*(1-(r0/a)^4)-(r0/a)^2*(1+(1+v)*log(a/r0)));
    F2 = 1/4 * (1 - (b./r).^2.*(1+2*log(r/b)));
    F3 = b./(4*r) .* ( ( (b./r).^2 + 1 ).*log(r/b) + (b./r).^2 - 1);
    F5 = 1/2 * (1 - (b./r).^2);
    F6 = b./(4*r) .* ( (b./r).^2 - 1 + 2*log(r/b) );
    F8 = 1/2 * (1+v+(1-v)*(b./r).^2);
    F9 = b./r .* (1/2*(1+v)*(log(r/b)) + 1/4*(1-v)*(1-(b./r).^2));
    

    bracket = zeros(size(r));
    bracket(r > r0) = 1;

    ratio = r0./r;
    G11 = 1/64 * (1 + 4*ratio.^2 - 5*ratio.^4 - 4*ratio.^2.*(2+ratio.^2).*log(1./ratio) ) .* bracket;
    G14 = 1/16 * (1 - (r0./r).^4  - 4*(r0./r).^2 .* log(r/r0)) .* bracket;
    G17 = 1/4 * (1 - ((1-v)/4)*(1-ratio.^4) - (ratio).^2.*(1+(1+v)*log(1./ratio))) .* bracket;

    % Mrb: reaction moment per unit length on inner edge
    Mrb = -q*a^2/C8 * (C9*(a^2-r0^2)/(2*a*b) - L17);
    % Qb: shear reaction force per unit length on inner edge
    Qb = q/2/b * (a^2 - r0^2);

    % case 2 general deflection equation p463, with first two terms
    % ommitted since yb and theta_b are zero for case 2L
    y_times_D = Mrb * r.^2 .* F2 ...
              + Qb * r.^3 .* F3 ...
              - q  * r.^4 .* G11;
    
    theta_times_D = Mrb * r .* F5 ...
                    + Qb * r.^2 .* F6 ...
                    - q * r.^3 .* G14;

    Mr = Mrb * F8 ...
        + Qb * r .* F9 ...
        - q * r.^2 .* G17;
    Mt = theta_times_D * (1-nu^2) ./ r ...
        + nu * Mr;
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
    w_nondim = y_times_D * 2*pi/(P*a^2);
    Mr_nondim = Mr * 2*pi/P;
    Mt_nondim = Mt * 2*pi/P;
    Q_nondim = Q * 2*pi/(P*a);

    % use sign convention that positive F_heave results in positive
    % deflection at outer end and positive moment at inner end
    w_nondim = -w_nondim;
    Mr_nondim = -Mr_nondim;
    Mt_nondim = -Mt_nondim;
    Q_nondim = -Q_nondim;
end