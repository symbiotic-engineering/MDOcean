
%syms a1 n
%get_C(0, n, a1)
%generate_BC_asdf(0,0,1,1,0,n,a1,4)

function BC = generate_BC(LHS_d_region, LHS_phi_region, RHS_d_region, RHS_phi_region, both_Z_region, m, a, N)

    % pass in five binary bits here for which region (outer/inner) is used
    % for which pieces of the matching condition
    % false = 0 = out, true = 1 = in
    % expected correct answer: 0, 0, 1, 1, 0
    
    % LHS (left hand side)  is defined as the side where orthogonality happens
    % RHS (right hand side) is defined as the side where coupling integral happens
    
    Z_m(m) = get_Z(both_Z_region, m, a);
    
    syms z h
    % LHS
    integrand_LHS = get_dphi_p_dr(LHS_phi_region, a) * Z_m;
    C_LHS(m) = get_C(LHS_phi_region, m, a);
    dR_dr_LHS(m) = get_dR_dr(LHS_phi_region, m, a);
    LHS = int(integrand_LHS, z, -h, LHS_d_region) + ...
        (h - LHS_d_region) * C_LHS * dR_dr_LHS.';
    
    % RHS
    j = sym('j'); % summation index
    integrand_RHS = get_dphi_p_dr(RHS_phi_region, a) * Z_m;
    Z_n(j) = get_Z(RHS_phi_region, j, a);
    coupling_integral(j,m) = int(Z_n(j) * Z_m(m), z, -h, RHS_d_region);
    C_RHS(j) = get_C(RHS_phi_region, j, a);
    dR_dr_RHS(j) = get_dR_dr(RHS_phi_region, j, a);
    sum_RHS(j,m) = C_RHS(j) * dR_dr_RHS(j).' * coupling_integral(j,m);
    RHS = int(integrand_RHS, z, -h, RHS_d_region) + symsum(sum_RHS, j, 0, N);
    
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

function lambda = lambda(j,region)
    syms d1 d2 h
    % eq 4
    if strcmp(region,'i1')
        di = d1;
    else
        di = d2;
    end
    lambda(j) = j*pi/(h-di);
end

function Z = get_Z(in_out, j, radius)
    syms z h
    region = get_region_from_in_out_radius(in_out,radius);

    if strcmp(region,'i1') || strcmp(region,'i2')
        % equation 9
        Z = piecewise(j==0, 1, ...
                      j>=1, sqrt(2)*cos(lambda(j,region)*(z+h)) ...
                      );
    elseif strcmp(region,'e')
        syms m0 m_k(dv)
        % dv means dummy variable: due to weird symbolic syntax it's necessary
        % to create the symfun as a function of the dummy and then replace it with j, 
        % instead of directly creating it as a function of j
        m_k(j) = subs(m_k,dv,j);

        % eq 2.34 in analytical methods book, also eq 16 in Seah and Yeung 2006 
        N_k(j) = piecewise(j==0, 1/2*(1+sinh(2*m0*h)/(2*m0*h)), ...
                           j>=1, 1/2*(1+sin(2*m_k(j)*h)/(2*m_k(j)*h)) ...
                           );
        % equation 14
        Z(j) = piecewise(j==0, 1/sqrt(N_k(j)) * cosh(m0 * (z+h)), ...
                         j>=1, 1/sqrt(N_k(j)) * cos(m_k(j) * (z+h)) ...
                         );
    end
end

function C = get_C(in_out, j, radius)
    region = get_region_from_in_out_radius(in_out,radius);

    % dv means dummy variable: due to weird symbolic syntax it's necessary
    % to create the symfun as a function of the dummy and then replace it with j, 
    % instead of directly creating it as a function of j
    syms C_1n_1(dv) C_1n_2(dv) C_2n_1(dv) B_k(dv)
    C_1n_1(j) = subs(C_1n_1,dv,j);
    C_1n_2(j) = subs(C_1n_2,dv,j);
    C_2n_1(j) = subs(C_2n_1,dv,j); 
    B_k(j)    = subs(B_k,   dv,j);

    if strcmp(region,'i1')
        C(j) = [C_1n_1(j), C_2n_1(j)];
    elseif strcmp(region,'i2')
        C(j) = [C_1n_2(j), 0];
    elseif strcmp(region,'e')
        C(j) = B_k(j);
    end

end

function dR_dr = get_dR_dr(in_out, j, radius)
    region = get_region_from_in_out_radius(in_out,radius);

    syms r
    R = get_R(region, j);
    dR_dr = subs( diff(R,r), r, radius);

end

function R = get_R(region, j)
    syms r a2 m0
   
    % dv means dummy variable: due to weird symbolic syntax it's necessary
    % to create the symfun as a function of the dummy and then replace it with j, 
    % instead of directly creating it as a function of j
    syms R(dv) m_k(dv)
    R(j) = subs(R,dv,j);
    m_k(j) = subs(m_k,dv,j);

    if strcmp(region,'i1') || strcmp(region,'i2')
        % eq 7
        R_1n(j) = piecewise(j==0, 1/2, j>=1, ...
                            besseli(0,lambda(j,region)*r)/besseli(0,lambda(j,region)*a2));
        if strcmp(region,'i1')
            % eq 8
            R_2n(j) = sym(0);
        else
            R_2n(j) = piecewise(j==0, 1/2*log(r/a2), ...
                                j>=1, besselk(0,lambda(j,region)*r)/besselk(0,lambda(j,region)*a2));
        end
        R = [R_1n, R_2n];

    elseif strcmp(region,'e')
        % eq 13
        Lambda_k(j) = piecewise(j==0, besselh(0,1,m0*r)/besselh(0,1,m0*a2), ...
                                j>=1, besselk(0,m_k(j)*r)/besselk(0,m_k(j)*a2));
        R = Lambda_k;
    end
end

function region = get_region_from_in_out_radius(in_out,radius)
    syms a1 a2
    if radius==a1 && in_out==1
        region = 'i1';
    elseif (radius==a1 && in_out==0) || (radius==a2 && in_out==1)
        region = 'i2';
    elseif radius==a2 && in_out==0
        region = 'e';
    else
        error('invalid radius or in/out value')
    end
end
