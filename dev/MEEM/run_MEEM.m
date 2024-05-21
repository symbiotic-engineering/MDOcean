function [mu_nondim, lambda_nondim] = run_MEEM(heaving_IC, heaving_OC, auto_BCs, ...
                                               N_num, M_num, K_num, ...
                                               a1_mat, a2_mat, d1_mat, ...
                                               d2_mat, h_mat, m0_mat, ...
                                               spatial_res, show_A, plot_phi)

    % Runs MEEM with potentially many numeric inputs (geometries and/or freqencies)
    % while keeping the same setup (number of harmonics, which cylinder is heaving).
    % Each numeric matrix is iterated through at the same time. This code
    % does not find the combinations of all possible inputs. If you want this, 
    % pre-generate the combinations (ie with meshgrid) and pass the meshes in here).

    % check valid geometry
    assert(all( d1_mat < h_mat,  'all') && ...
           all( d2_mat < h_mat,  'all') && ...
           all( a1_mat < a2_mat, 'all') && ...
           all( d2_mat < d1_mat, 'all'), 'Invalid geometry')

    % check input dimensions
    numel_numeric_inputs = [numel(a1_mat), numel(a2_mat), numel(d1_mat), numel(d2_mat), ...
                     numel(h_mat), numel(m0_mat)];
    numel_nonscalar = numel_numeric_inputs(numel_numeric_inputs ~= 1);

    if isempty(numel_nonscalar) % all scalars
        num_runs = 1;
    else                        % some vectors/matrices
        num_runs = unique(numel_nonscalar);
        assert( isscalar(num_runs), ...
            'All nonscalar numeric inputs must have same number of elements');
    end

    % generate filename based on setup
    if heaving_IC && heaving_OC
        heaving = 'both';
    elseif heaving_IC
        heaving = 'inner';
    else
        heaving = 'outer';
    end
    fname = ['N' num2str(N_num) '_M' num2str(M_num) '_K' num2str(K_num) '_heaving_' heaving];

    % generate matlab functions from symbolic equations if needed
    if ~exist(['dev/MEEM/generated/A_b_matrix_' fname],'file') || ...
       ~exist(['dev/MEEM/generated/hydro_potential_velocity_fields_' fname],'file')
        create_symbolic_expressions(heaving_IC, heaving_OC, auto_BCs, N_num, M_num, K_num, fname)
    end  

    % for loop for numeric inputs
    mu_nondim = zeros(1,num_runs);
    lambda_nondim = zeros(1,num_runs);
    
    for i = 1:num_runs
        % index into matrix inputs with i and into scalar inputs with 1
        idxs = numel_numeric_inputs;
        idxs(numel_numeric_inputs~=1) = i;

        a1_num = a1_mat(idxs(1));
        a2_num = a2_mat(idxs(2));
        d1_num = d1_mat(idxs(3));
        d2_num = d2_mat(idxs(4));
        h_num  =  h_mat(idxs(5));
        m0_num = m0_mat(idxs(6));

        [mu_nondim(i), lambda_nondim(i)] = compute_and_plot(a1_num, a2_num, d1_num, d2_num, ...
                                                                h_num, m0_num, spatial_res, ...
                                                                K_num, show_A, plot_phi, fname);
    end
end

function [mu_nondim, lambda_nondim] = compute_and_plot(a1_num, a2_num, d1_num, d2_num, h_num, m0_num, spatial_res, K_num, show_A, plot_phi, fname)
    % solve for m_k from m_0 and h    
    
    [x_cell, m_k_cell] = compute_eigen_coeffs(a1_num,a2_num,d1_num,d2_num,h_num,m0_num,K_num,show_A,fname);

    hydro_fname = ['hydro_potential_velocity_fields_' fname];
    if plot_phi
        % generate R Z mesh
        r_vec = linspace(2*a2_num/spatial_res,2*a2_num,spatial_res);
        a_eps = 1e-4;
        r_vec = sort([r_vec,a1_num*(1-a_eps), a1_num*(1+a_eps), ...
                            a2_num*(1-a_eps), a2_num*(1+a_eps)]);
        z_vec = linspace(-h_num,0,spatial_res);
        [R,Z] = meshgrid(r_vec,z_vec);
        
        % get hydro coeffs and potential velocity fields
        [hydro_nondim_num,phi_i1_num,phi_i2_num,...
        phi_e_num,phi_p_i1_num,phi_p_i2_num,...
        phi_h_i1_num,phi_h_i2_num,v_1_r_num,...
        v_1_z_num,v_2_r_num,v_2_z_num,v_e_r_num,...
        v_e_z_num] = feval(hydro_fname, a1_num,a2_num,d1_num,d2_num,h_num,...
                                                    m0_num,m_k_cell{:},x_cell{:},R,Z);
        [phi_p_i1_num,phi_p_i2_num] = fix_scalars(size(R),phi_p_i1_num,phi_p_i2_num);
        assemble_plot_pot_vel_fields(a1_num,a2_num,d1_num,d2_num,R,Z,phi_i1_num,phi_i2_num,phi_e_num,...
                                     phi_p_i1_num,phi_p_i2_num,phi_h_i1_num,phi_h_i2_num,...
                                     v_1_r_num,v_1_z_num,v_2_r_num,v_2_z_num,v_e_r_num,v_e_z_num,fname);
    else
        % just get hydro coeffs
        hydro_nondim_num = feval(hydro_fname,a1_num,a2_num,d1_num,d2_num,h_num,...
                                                    m0_num,m_k_cell{:},x_cell{:},0,0);
    end

    g = 9.806;
    omega = sqrt(g * m0_num * tanh(m0_num*h_num));
    mu_nondim = real(hydro_nondim_num);
    lambda_nondim = imag(hydro_nondim_num) / omega;

end

function [varargout] = fix_scalars(desired_size, varargin)
% ensures that all arguments passed have desired_size. If they are scalar,
% converts them to the desired size by multiplying by appropriate one matrix.

    varargout = cell(1,numel(varargin));
    for i = 1:numel(varargin)
        varargout{i} = ones(desired_size) .* varargin{i};
    end
end

function [x_cell, m_k_cell] = compute_eigen_coeffs(a1_num,a2_num,d1_num,d2_num,h_num,m0_num,K_num,show_A,fname)

    m_k_num = zeros(1,K_num);
    % using tand instead of tan because finite precision of pi means
    % inconsistent behaviour around tan(pi/2)
    eqn = @(m_k_h_deg) (m_k_h_deg * pi/180) .* tand(m_k_h_deg) + m0_num * h_num * tanh(m0_num * h_num);

    for k_num = 1:K_num
        bounds = 180 * [k_num-1/2, k_num];
        % apply tweak determined experimentally to be the smallest tweak
        % that doesn't return +-Inf - see p17 of notebook
        expo = floor( log(bounds(1)) / log(2) );
        bound_tweak = eps * (2^(expo - 1) + 1);
        bounds(1) = bounds(1) + bound_tweak;

        m_k_h_deg = fzero(eqn, bounds);
        m_k_num(k_num) = m_k_h_deg * pi/180 / h_num;
    end    

    m_k_cell = num2cell(m_k_num);

    % get A and b matrices
    Ab_fname = ['A_b_matrix_' fname];
    [A_num, b_num] = feval(Ab_fname, a1_num, a2_num, d1_num, d2_num,...
                            h_num, m0_num, m_k_cell{:});
    
    % show A matrix values
    if show_A
        cond(A_num)
        signed_log(real(A_num))
        title('Real(A)')
        
        signed_log(imag(A_num))
        title('Imag(A)')
    end

    % solve for x
    x = full(A_num\b_num);
    x_cell = num2cell(x);
end

function create_symbolic_expressions(heaving_IC, heaving_OC, auto_BCs, N_num, M_num, K_num, fname)

    syms r z real                                           % coordinates
    syms R_1n_1(n) R_1m_2(m) R_2m_2(m) Z_n_i1(n)  ...       % basis functions
            Z_m_i2(m) Lambda_k(k) N_k(n) Z_k_e(k) m_k(k)
    syms C_1n_1(n) C_1m_2(m) C_2n_1(n) C_2m_2(m) B_k(k)     % unknown coefficients              
    syms h a1 a2 d1 d2 m0 real positive                     % constants
    syms n m k real integer                                 % counter index
    assumeAlso(n >= 0)
    assumeAlso(m >= 0)
    assumeAlso(k >= 0)
    
    % equation numbers refer to Chau & Yeung 2010 unless otherwise noted
    
    %% setup analytical boundary value problem equations
    % eq 4
    lambda_n1(n) = n*pi/(h-d1);
    lambda_m2(m) = m*pi/(h-d2);
    
    % eq 5
    if heaving_IC
        phi_p_i1 = 1/(2*(h-d1)) * ((z+h)^2 - r^2/2);
    else
        phi_p_i1 = 0;
    end
    if heaving_OC
        phi_p_i2 = 1/(2*(h-d2)) * ((z+h)^2 - r^2/2);
    else
        phi_p_i2 = 0;
    end
    
    % eq 7
    R_1n_1(n) = piecewise(n==0, 1/2, n>=1, ...
                        besseli(0,lambda_n1(n)*r)/besseli(0,lambda_n1(n)*a2));
    R_1m_2(m) = piecewise(m==0, 1/2, m>=1, ...
                        besseli(0,lambda_m2(m)*r)/besseli(0,lambda_m2(m)*a2));
    
    % eq 8
    R_2n_1(n) = sym(0);
    R_2m_2(m) = piecewise(m==0, 1/2*log(r/a2), ...
        m>=1, besselk(0,lambda_m2(m)*r)/besselk(0,lambda_m2(m)*a2));
    
    % eq 9
    Z_n_i1(n) = piecewise(n==0, 1, n>=1, sqrt(2)*cos(lambda_n1(n)*(z+h)));
    Z_m_i2(m) = piecewise(m==0, 1, m>=1, sqrt(2)*cos(lambda_m2(m)*(z+h)));
    
    % eq 13
    Lambda_k(k) = piecewise(k==0, besselh(0,1,m0*r)/besselh(0,1,m0*a2), ...
        k>=1, besselk(0,m_k(k)*r)/besselk(0,m_k(k)*a2));
    
    % eq 2.34 in analytical methods book, also eq 16 in Seah and Yeung 2006 
    N_k(k) = piecewise(k==0, 1/2*(1+sinh(2*m0*h)/(2*m0*h)), ...
                       k>=1, 1/2*(1+sin(2*m_k(k)*h)/(2*m_k(k)*h)) );
    
    % eq 14
    Z_k_e(k) = piecewise(k==0, 1/sqrt(N_k(k)) * cosh(m0 * (z+h)), ...
                         k>=1, 1/sqrt(N_k(k)) * cos(m_k(k) * (z+h)));
    
    
    % coupling integrals
    C_nm(n,m) = simplify(int(Z_m_i2(m) * Z_n_i1(n), z, -h, -d1));
    C_mk(m,k) = simplify(int(Z_m_i2(m) * Z_k_e(k),  z, -h, -d2));
    
    
    dz_1 = h - d1;
    dz_2 = h - d2;
    
    % potential matching
    % equation 22 in old 1981 paper, applied to boundary 2-e
    match_2e_potential = dz_2 * ( C_1m_2(m)*subs(R_1m_2(m),r,a2) ...
        + C_2m_2(m)*subs(R_2m_2(m),r,a2) ) + int( subs(phi_p_i2,r,a2) * Z_m_i2(m), z, -h, -d2) == ...
        symsum( B_k(k) * subs(Lambda_k,r,a2) * C_mk(m,k), k, 0, K_num);
    
    % equation 22 in old 1981 paper, applied to boundary 1-2
    match_12_potential = C_1n_1(n) * subs(R_1n_1(n),r,a1) * dz_1 == ...
        symsum(( (C_1m_2(m) * subs(R_1m_2(m),r,a1) + ...
          C_2m_2(m) * subs(R_2m_2(m),r,a1) ) * C_nm(n,m) ),m,0,M_num) + ...
        int(subs(phi_p_i2,r,a1)*Z_n_i1(n) - subs(phi_p_i1,r,a1)*Z_n_i1(n), z, -h, -d1);
    
    % velocity matching
    % equation 23 in old 1981 paper, applied to boundary 2-e
    match_2e_velocity = B_k(k) * subs(diff(Lambda_k(k), r), r, a2) * (h) == ...
        symsum((C_1m_2(m) * subs(diff(R_1m_2(m), r), r, a2) + ...
         C_2m_2(m) * subs(diff(R_2m_2(m), r), r, a2) ) * C_mk(m,k), m, 0, M_num) + ...
        int(subs(diff(phi_p_i2,r),r,a2) * Z_k_e(k), z, -h, -d2 );
    
    % equation 23 in old 1981 paper, applied to boundary 1-2
    %     use large region (-h to -d2) and multiply by Z_m_i2
    match_12_velocity = (h-d2) * ( C_1m_2(m) * subs(diff(R_1m_2(m), r), r, a1) + ...
     C_2m_2(m) * subs(diff(R_2m_2(m), r), r, a1) ) + int( subs(diff(phi_p_i2,r),r,a1) * Z_m_i2(m), z, -h, -d2 ) == ...
     symsum((C_1n_1(n) * subs(diff(R_1n_1(n), r), r, a1) + ...
     C_2n_1(n) * subs(diff(R_2n_1(n), r), r, a1) ) * C_nm(n,m), n, 0, N_num) + ...
     int( subs(diff(phi_p_i1,r),r,a1) * Z_m_i2(m), z, -h, -d1 );
    
    if auto_BCs
        error('auto_BCs are broken right now, please set auto_BCs to false')
        % enter a number 0 to 31 for the combination
        % expected best is 6 = 0 0 1 1 0
        num = 6;
        % done: 29,25,23,22,21,19,18,17,13,9,7,6,5,3,2,1
        % all Nan: 30,28,27,26,24,20,16,15,14,12,11,10,8,4,0
        binary = dec2bin(num, 5);
        
        % Assign each bit to a separate variable
        bit1 = str2double(binary(1));
        bit2 = str2double(binary(2));
        bit3 = str2double(binary(3));
        bit4 = str2double(binary(4));
        bit5 = str2double(binary(5));
        
        match_12_velocity = generate_BC(bit1,bit2,bit3,bit4,bit5,m,a1,N_num);
        match_2e_velocity = generate_BC(bit1,bit2,bit3,bit4,bit5,k,a2,M_num);
    end
    
    % get potential field
    % eq 6
    phi_h_n_i1(n) = (C_1n_1(n) * R_1n_1(n) + C_2n_1(n) * R_2n_1(n)) * Z_n_i1(n);
    phi_h_m_i2(n) = (C_1m_2(m) * R_1m_2(m) + C_2m_2(m) * R_2m_2(m)) * Z_m_i2(m);
    
    % eq 12
    phi_e_k(k) = B_k(k) .* Lambda_k(k) .* Z_k_e(k);
    
    % simplify
    phi_h_n_i1(n) = simplify(phi_h_n_i1(n),'IgnoreAnalyticConstraints',true);
    phi_h_m_i2(m) = simplify(phi_h_m_i2(m),'IgnoreAnalyticConstraints',true);
    phi_e_k(k)    = simplify(phi_e_k(k),   'IgnoreAnalyticConstraints',true);
    
    % sum all N terms
    syms N M K real positive integer
    phi_h_i1 = symsum(phi_h_n_i1, n, 0, N);
    phi_h_i2 = symsum(phi_h_m_i2, m, 0, M);
    phi_i1 = phi_p_i1 + phi_h_i1;
    phi_i2 = phi_p_i2 + phi_h_i2;
    phi_e  = symsum(phi_e_k, k, 0, K);
    
    % derivatives to geometric variables
    % hndInfDepth = limit(hnd,h,inf)
    % geom_vars = [a1,a2,d1,d2,h];
    % hydro_terms = [mu,lambda];
    % J = jacobian(hydro_terms,geom_vars);
    % matlabFunction(J,'File','generated/derivatives','Vars',symvar(J))
    % J = J(n) % turn symfun into sym so I can index matrix elements
    % spy(J)
    % xticklabels({'','a_1','a_2','d_1','d_2','h'})
    % yticklabels({'','\mu','\lambda'})
    
    % get symbolic velocity
    v_1r = diff(phi_i1,r);
    v_1z = diff(phi_i1,z);
    v_2r = diff(phi_i2,r);
    v_2z = diff(phi_i2,z);
    v_er = diff(phi_e,r);
    v_ez = diff(phi_e,z);
    
    % get symbolic hydro coeffs
    % equation 28 in old 1981 paper
    % OC = heaving outer cylinder
    % IC = heaving inner cylinder
    integrand_OC = subs(r * phi_i2 * v_2z, z, -d2);
    integrand_IC = subs(r * phi_i1 * v_1z, z, -d1);
    hydro_OC = h^3 * int(integrand_OC,r,a1,a2);
    hydro_IC = h^3 * int(integrand_IC,r,0,a1);
    
    if heaving_OC && heaving_IC
        hydro_over_rho = hydro_OC + hydro_IC;
        a_normalize = a2;
    elseif heaving_OC
        hydro_over_rho = hydro_OC;
        a_normalize = a2;
    elseif heaving_IC
        hydro_over_rho = hydro_IC;
        a_normalize = a1;
    else
        error('neither body is setup to move')
    end
    hydro_nondim = hydro_over_rho / (pi * a_normalize^3);
    
    % The above was general for any N,M,K. Now we sub in particular N_num, M_num, K_num.
    % These substituted variables will be called _NMK.
    
    N_vec = (0:N_num);
    M_vec = (0:M_num);
    K_vec = (0:K_num);
    
    % set up equations
    eqns = [subs(match_12_potential,n,N_vec), subs(match_2e_potential,m,M_vec),...
            subs(match_12_velocity, m,M_vec), subs(match_2e_velocity,k,K_vec) ].';
    
    unknowns = [C_1n_1(N_vec) C_1m_2(M_vec) C_2m_2(M_vec) B_k(K_vec)];
    syms C_1n_1_const [1 N_num+1]
    syms C_1m_2_const C_2m_2_const [1 M_num+1]
    syms B_k_const [1 K_num+1]
    unknowns_const = [C_1n_1_const C_1m_2_const C_2m_2_const B_k_const];
    syms m_k_const [1 K_num] real positive
    
    var = {phi_i1,phi_i2,phi_e,phi_p_i1,phi_p_i2,phi_h_i1,phi_h_i2,...
        v_1r,v_1z,v_2r,v_2z,v_er,v_ez,hydro_nondim,eqns};
    for i = 1:length(var)
        var{i} = subs(var{i}, [N M K], [N_num M_num K_num]);
        var{i} = subs(var{i}, [unknowns m_k(1:K_num)], [unknowns_const m_k_const]);
    end
    
    phi_i1_NMK = var{1};
    phi_i2_NMK = var{2};
    phi_e_NMK = var{3};
    phi_p_i1_NMK = var{4};
    phi_p_i2_NMK = var{5};
    phi_h_i1_NMK = var{6};
    phi_h_i2_NMK = var{7};
    v_1_r_NMK = var{8};
    v_1_z_NMK = var{9};
    v_2_r_NMK = var{10};
    v_2_z_NMK = var{11};
    v_e_r_NMK = var{12};
    v_e_z_NMK = var{13};
    hydro_nondim_NMK = var{14};
    eqns_NMK = var{15};
    
    % generate linear system
    [A,b] = equationsToMatrix(eqns_NMK,unknowns_const);
    
    % visualize sparsity pattern of A
    clf
    spy(A)
    hold on
    widths = [N_num+1, M_num+1, M_num+1, K_num+1];
    bars = .5 + [widths(1), sum(widths(1:2)), sum(widths(1:3))];
    full_line = [0,1+sum(widths)];
    
    plot(bars(1)*[1,1], full_line, 'k') % first vertical
    plot(bars(2)*[1,1], full_line, 'k') % second vertical
    plot(bars(3)*[1,1], full_line, 'k') % third vertical
    
    plot(full_line, bars(1)*[1,1], 'k') % first horizontal
    plot(full_line, bars(2)*[1,1], 'k') % first horizontal
    plot(full_line, bars(3)*[1,1], 'k') % first horizontal
    
    matlabFunction(A,b,'File',['dev/MEEM/generated/A_b_matrix_' fname], 'Vars',[a1,a2,d1,d2,h,m0,m_k_const]);
    
    matlabFunction(hydro_nondim_NMK,phi_i1_NMK,phi_i2_NMK,phi_e_NMK,phi_p_i1_NMK,phi_p_i2_NMK,...
                    phi_h_i1_NMK,phi_h_i2_NMK,v_1_r_NMK,v_1_z_NMK,v_2_r_NMK,v_2_z_NMK,v_e_r_NMK,v_e_z_NMK,...
                    'File',['dev/MEEM/generated/hydro_potential_velocity_fields_' fname],...
                    'Vars',[a1,a2,d1,d2,h,m0,m_k_const,unknowns_const,r,z]);

end

function assemble_plot_pot_vel_fields(a1_num,a2_num,d1_num,d2_num,R,Z,...
                                     phi_i1_num,phi_i2_num,phi_e_num,...
                                     phi_p_i1_num,phi_p_i2_num,phi_h_i1_num,phi_h_i2_num,...
                                     v_1_r_num,v_1_z_num,v_2_r_num,v_2_z_num,v_e_r_num,v_e_z_num,...
                                     fname)
    % assemble total phi based on phi in each region
    regione = R > a2_num;
    region1 = R <= a1_num & Z < -d1_num;
    region2 = R > a1_num & R <= a2_num & Z < -d2_num;
    
    phi = NaN(size(R));
    phi(region1) = phi_i1_num(region1);
    phi(region2) = phi_i2_num(region2);
    phi(regione) = phi_e_num(regione);
    
    phiH = NaN(size(R));
    phiH(region1) = phi_h_i1_num(region1);
    phiH(region2) = phi_h_i2_num(region2);
    
    phiP = NaN(size(R));
    phiP(region1) = phi_p_i1_num(region1);
    phiP(region2) = phi_p_i2_num(region2);
    
    v_r = NaN(size(R));
    v_r(region1) = v_1_r_num(region1);
    v_r(region2) = v_2_r_num(region2);
    v_r(regione) = v_e_r_num(regione);
    
    v_z = NaN(size(R));
    v_z(region1) = v_1_z_num(region1);
    v_z(region2) = v_2_z_num(region2);
    v_z(regione) = v_e_z_num(regione);

    region_body = ~region1 & ~region2 & ~regione;
    plot_potential(phi,R,Z,region_body,'Total Potential');

    date_string = [char(datetime('now','Format','yyyy-MM-dd__h-mma')) '_'];
    folder = ['dev' filesep 'MEEM' filesep 'MEEM_figs' filesep date_string fname];
    savefig(folder)
    plot_potential(phiH,R,Z,region_body,'Homogeneous Potential');
    plot_potential(phiP,R,Z,region_body,'Particular Potential');
    plot_potential(v_r,R,Z,region_body,'Radial Velocity')
    plot_potential(v_z,R,Z,region_body,'Vertical Velocity')
    plot_velocity(real(v_r),real(v_z),R,Z);

    z_cutoffs_potential = [-d1_num -d1_num -d2_num -d2_num];
    z_cutoffs_velocity = [-d1_num -d2_num -d2_num 0];

    plot_matching(phi_i1_num,phi_i2_num,phi_e_num,a1_num,a2_num,R,Z,'potential',z_cutoffs_potential)

    plot_matching(v_1_r_num,v_2_r_num,v_e_r_num,a1_num,a2_num,R,Z,'Radial Velocity',z_cutoffs_velocity)
    plot([-d1_num -d2_num],[0 0],'m','DisplayName','No flux BC at a1')
    plot([-d2_num 0]      ,[0 0],'c--','DisplayName','No flux BC at a2')
    savefig([folder '_VelMatch'])
end

function [] = plot_potential(phi,R,Z,region_body,name)
    figure
    minphi = min(real(phi),[],'all');
    maxphi = max(real(phi),[],'all');

    num_levels = 20;
    levels = linspace(minphi,maxphi,num_levels);

    subplot 121
    [c,h_fig] = contourf(R,Z,real(phi),levels);
    clabel(c,h_fig)
    xlabel('R')
    ylabel('Z')
    title([name ' - Real'])
    colorbar;
    
    imag_phi = imag(phi);
    if numel(unique(imag_phi)) > 1
        imag_phi(region_body) = NaN;
        subplot 122
        [c,h_fig] = contourf(R,Z,imag_phi,num_levels);
        clabel(c,h_fig)
        xlabel('R')
        ylabel('Z')
        title([name ' - Imaginary'])
        colorbar
    end
end

function plot_matching(phi1,phi2,phie,a1,a2,R,Z,name,z_cutoffs)
    % at R = a1
    idx_a1_minus = R < a1 & R == max(R(R < a1)) & Z < z_cutoffs(1);
    idx_a1_plus  = R > a1 & R == min(R(R > a1)) & Z < z_cutoffs(2);
    phi1_a1 = abs(phi1(idx_a1_minus));
    phi2_a1 = abs(phi2(idx_a1_plus));

    % at R = a2
    idx_a2_minus = R < a2 & R == max(R(R < a2)) & Z < z_cutoffs(3);
    idx_a2_plus  = R > a2 & R == min(R(R > a2)) & Z < z_cutoffs(4);
    phi2_a2 = abs(phi2(idx_a2_minus));
    phie_a2 = abs(phie(idx_a2_plus));

    % plot
    figure
    plot(Z(idx_a1_minus),phi1_a1,'r--', ...
         Z(idx_a1_plus), phi2_a1,'m-')
    hold on
    plot(Z(idx_a2_minus),phi2_a2,'b-', ...
         Z(idx_a2_plus), phie_a2,'c--')
    legend([name '_1 at a_1'],[name '_2 at a_1'],[name '_2 at a_2'],[name '_e at a_2'])
    xlabel('Z')
    ylabel(['|' name '|'])
    title([name ' Matching'])
end

function plot_velocity(v_r,v_z,R,Z)
    v_tot = sqrt(v_r.^2 + v_z.^2);

    num_levels = 10;
    levels = linspace(1,6,num_levels);

    figure
    [c,h_fig] = contourf(R,Z,v_tot,levels);
    clabel(c,h_fig)
    xlabel('R')
    ylabel('Z')
    title('Velocity')
    colorbar
    hold on
    quiver(R,Z,v_r./v_tot,v_z./v_tot)
end
