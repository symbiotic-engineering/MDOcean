function [w_nondim,Mr_nondim,abcd] = concentrated_plate_nondim(lam,nu,theta,rho,N)
% Boedo and Prantil 1998: corrected solution of clamped ring plate with edge point load
% https://ascelibrary.org/doi/epdf/10.1061/%28ASCE%290733-9399%281998%29124%3A6%28696%29

    % lam: aspect ratio: b/a = inner radius/outer radius
    % nu: poisson ratio
    % theta: angular coordinate
    % rho: nondim radial coordinate
    % N: number of terms to compute in infinite sum

    assert(isscalar(lam))
    assert(isscalar(nu))
    assert(isscalar(N))
    % vector inputs only allowed for theta and rho

    [RHO,n] = meshgrid(rho,2:N);
    d0 = -1/4;
    c0_num = (1 + 2*log(lam)) * (1+nu) - (3+nu);
    c0_den = lam^(-1)*(1+nu) + lam*(1-nu);
    c0 = lam/4 * c0_num / c0_den;
    b0 = c0/2 * (1-nu)/(1+nu) + (3+nu)/(8*(1+nu));
    a0 = -lam^2*b0 - c0*log(lam) - d0*lam^2*log(lam);
    
    d1 = 1/2;
    b1 = -1/4 * (1+nu+lam^2*(1-nu)) / (3+nu+lam^4*(1-nu));
    c1 = -b1*(3+nu)/(1-nu) - 1/4 * (1+nu)/(1-nu);
    a1 = -b1*lam^2 - c1*lam^(-2) - d1*log(lam);
    
    A = (3+nu)/(1-nu);
    B = (1-lam^2)^2 * (n.^2-1) + (lam.^(-2*n+2)+A) .* (lam.^(2*n+2)+A);
    
    dn_num = (1-lam^2)*(n-1) + lam.^( 2*n+2) + A;
    bn_num = (1-lam^2)*(n+1) - lam.^(-2*n+2) - A;
    dn_denom = B.*n.*(n-1)*(1-nu);
    bn_denom = B.*n.*(n+1)*(1-nu);
    dn = -dn_num ./ dn_denom;
    bn =  bn_num ./ bn_denom;
    an = -lam^2 * (bn .* (n+1)./n + dn .* lam.^(-2*n)./n);
    cn = -lam^2 * (dn .* (n-1)./n - bn .* lam.^( 2*n)./n);
    
    if length(rho)==1
        abcd  = [a0 b0 c0 d0;
                 a1 b1 c1 d1;
                 an bn cn dn];
    else
        abcd = []; % don't bother outputting coeffs if they're all vectors
    end
    
    Rho_zero_terms = a0 + b0 * rho.^2 + c0 * log(rho) + d0 * rho.^2 .* log(rho);
    Rho_one_terms = a1*rho + b1 * rho.^3 + c1 * rho.^(-1) + d1 * rho .* log(rho);

    
    Rho_two_to_N_terms = an .* RHO.^n + bn .* RHO.^(n+2) + cn .* RHO.^(-n) + dn .* RHO.^(-n+2);
    
    w_nondim = Rho_zero_terms' + Rho_one_terms' * cos(theta) + Rho_two_to_N_terms' * cos(n(:,1)*theta);
    
    e0_d0_coeff = 3 + nu + 2*(1+nu)*log(rho);
    e0 = 2*b0*(1+nu) - c0*rho.^-2*(1-nu) + d0*(e0_d0_coeff);
    e1 = 2*b1*rho*(3+nu) + 2*c1*rho.^-3*(1-nu) + d1*rho.^-1*(1+nu);
    
    en_a_term = an .* RHO.^(n-2) .* n .* (n-1) * (1 - nu);
    en_b_term = bn .* RHO.^n .* (n+1) .* (n + 2 - nu*(n-2));
    en_c_term = cn .* RHO.^(-n-2) .* n .* (n+1) * (1-nu);
    en_d_term = dn .* RHO.^-n .* (n-1) .* (n-2 - nu*(n+2));
    en = en_a_term + en_b_term + en_c_term + en_d_term;
    Mr_nondim = e0' + e1'*cos(theta) + en'*cos(n(:,1)*theta);

    % use sign convention that positive input force results in negative
    % deflection at outer end and negative moment at inner end
    w_nondim = -w_nondim;
    Mr_nondim = -Mr_nondim;
end