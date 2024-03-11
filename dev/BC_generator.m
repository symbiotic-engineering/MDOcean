syms m a_1 a_2
N = 4;
BC_a1 = generate_BC(0,0,1,1,0,m,a_1,N)
BC_a2 = generate_BC(0,0,1,1,0,m,a_2,N)

function BC = generate_BC(LHS_d_region, LHS_phi_region, RHS_d_region, RHS_phi_region, both_Z_region, m, a, N)

    % pass in five binary bits here
    % false = 0 = out, true = 1 = in
    % expected correct answer: 0, 0, 1, 1, 0
    
    % LHS (left hand side)  is defined as the side where orthogonality happens
    % RHS (right hand side) is defined as the side where coupling integral happens
    
    Z_m = get_Z(both_Z_region, m, a);
    
    syms z h
    % LHS
    integrand_LHS = get_dphi_p_dr(LHS_phi_region, a) * Z_m;
    LHS = int(integrand_LHS, z, -h, LHS_d_region) + ...
        (h - LHS_d_region) * get_C(LHS_phi_region, m, a) .* get_dR_dr(LHS_phi_region, m, a);
    
    % RHS
    syms n % summation index
    integrand_RHS = get_dphi_p_dr(RHS_phi_region, a) * Z_m;
    Z_n = get_Z(RHS_phi_region, n, a);
    coupling_integral(n,m) = int(Z_n(n) * Z_m(m), z, -h, RHS_d_region);
    sum_RHS(n,m) = get_C(RHS_phi_region,n, a) .* get_dR_dr(RHS_phi_region,n, a) * coupling_integral(n,m);
    RHS = int(integrand_RHS, z, -h, RHS_d_region) + symsum(sum_RHS, n, 0, N);
    
    % equation
    BC = LHS == RHS;

end

function d_phi_p_dr = get_dphi_p_dr(in_out, radius)
    syms d1 d2 h z r
    region = get_region_from_in_out_radius(in_out,radius);
    if strcmp(region,'i1') || strcmp(region,'i2')
        if strcmp(region,'i1')
            di = d1;
        else
            di = d2;
        end
        % eq 5
        phi_p = 1/(2*(h-di)) * ((z+h)^2 - r^2/2);
        d_phi_p_dr = subs( diff(phi_p, r), r, radius);
    elseif strcmp(region,'e')
        d_phi_p_dr = 0;
    end

end

function lambda = lambda(n,region)
    syms d1 d2 h
    % eq 4
    if strcmp(region,'i1')
        di = d1;
    else
        di = d2;
    end
    lambda(n) = n*pi/(h-di);
end

function Z = get_Z(in_out, n, radius)
    syms z h
    region = get_region_from_in_out_radius(in_out,radius);

    if strcmp(region,'i1') || strcmp(region,'i2')
        % equation 9
        Z = piecewise(n==0, 1, ...
                      n>=1, sqrt(2)*cos(lambda(n,region)*(z+h)) ...
                      );
    elseif strcmp(region,'e')
        syms m0 m_k
        % eq 2.34 in analytical methods book, also eq 16 in Seah and Yeung 2006 
        N_k(n) = piecewise(n==0, 1/2*(1+sinh(2*m0*h)/(2*m0*h)), ...
                                  n>=1, 1/2*(1+sin(2*m_k(n)*h)/(2*m_k(n)*h)) ...
                                  );
        % equation 14
        Z = piecewise(n==0, 1/sqrt(N_k(n)) * cosh(m0 * (z+h)), ...
                      n>=1, 1/sqrt(N_k(n)) * cos(m_k(n) * (z+h)) ...
                      );
    end
end

function C = get_C(in_out, n, radius)
    region = get_region_from_in_out_radius(in_out,radius);
    syms C_1n_1 C_1n_2 C_2n_1 B_k
    symfun(C_1n_1,n); symfun(C_1n_2,n); symfun(C_2n_1,n); symfun(B_k,n);
    if strcmp(region,'i1')
        C = [C_1n_1, C_2n_1];
    elseif strcmp(region,'i2')
        C = [C_1n_2, 0];
    elseif strcmp(region,'e')
        C = B_k;
    end

end

function dR_dr = get_dR_dr(in_out, n, radius)
    region = get_region_from_in_out_radius(in_out,radius);

    syms r
    R = get_R(region, n);
    dR_dr = subs( diff(R,r), r, radius);

end

function R = get_R(region, n)
    syms r a2 m0 R m_k
    symfun(R,n); symfun(m_k,n);
    if strcmp(region,'i1') || strcmp(region,'i2')
        % eq 7
        R_1n(n) = piecewise(n==0, 1/2, n>=1, ...
                            besseli(0,lambda(n,region)*r)/besseli(0,lambda(n,region)*a2));
        if strcmp(region,'i1')
            % eq 8
            R_2n(n) = sym(0);
        else
            R_2n(n) = piecewise(n==0, 1/2*log(r/a2), ...
                                n>=1, besselk(0,lambda(n,region)*r)/besselk(0,lambda(n,region)*a2));
        end
        R = [R_1n, R_2n];

    elseif strcmp(region,'e')
        % eq 13
        Lambda_k(n) = piecewise(n==0, besselh(0,1,m0*r)/besselh(0,1,m0*a2), ...
                                n>=1, besselk(0,m_k(n)*r)/besselk(0,m_k(n)*a2));
        R = Lambda_k;
    end
end

function region = get_region_from_in_out_radius(in_out,radius)
    syms a_1 a_2
    if radius==a_1 && in_out==1
        region = 'i1';
    elseif (radius==a_1 && in_out==0) || (radius==a_2 && in_out==1)
        region = 'i2';
    elseif radius==a_2 && in_out==0
        region = 'e';
    else
        error('invalid radius or in/out value')
    end
end
