function [A, b, dAdp, dbdp] = lin_ineq_constraints(p, param_name)

MEEM = pi*p.harmonics / (p.besseli_argmax*2);

% A x < b

%    x1                      x2 x3                      x4 x5 x6 x7 x8 x9 x10 x11 x12 x13
A = [-p.D_d_over_D_s         0  0                        0  0  0  0  0  0  0   0   0   0; % g_lin_1
      p.D_f_in_over_D_s     -1  0                        0  0  0  0  0  0  0   0   0   0; % g_lin_2 
     -p.T_s_over_D_s         0  1                        0  0  0  0  0  0  0   0   0   0; % g_lin_3
      p.T_s_over_D_s,        0, 1/p.T_f_2_over_h_f - 1, -1  0  0  0  0  0  0   0   0   0; % g_lin_4 
      0                    MEEM 1                        0  0  0  0  0  0  0   0   0   0; % g_lin_5 
      MEEM + p.T_s_over_D_s, 0  0                        0  0  0  0  0  0  0   0   0   0; % g_lin_6 
     -p.h_d_over_D_s         0  0                        0  0  0  0  0  0 1e-3 0   0   0; % g_lin_7 
      0                      0  -.5/p.T_f_2_over_h_f,    0  0  0  0  0  0  0   1   0   0];% g_lin_8 

%     g1         g2    g3    g4   g5  g6   g7    g8
b = [-p.D_d_min -0.01 -0.01 -0.01 p.h p.h -0.01 -0.01]';


if nargout > 2 && nargin > 1
    dAdp = zeros(size(A));
    dbdp = zeros(size(b));
    if strcmp(param_name,'D_d_over_D_s')
        dAdp(1,1) = -1;
    elseif strcmp(param_name,'T_s_over_D_s')
        dAdp(3,1) = -1;
        dAdp(4,1) = 1;
        dAdp(6,1) = 1;
    elseif strcmp(param_name,'T_f_2_over_h_f')
        dAdp(4,3) = -1/p.T_f_2_over_h_f^2;
        dAdp(8,3) = .5/p.T_f_2_over_h_f^2;
    elseif strcmp(param_name,'harmonics')
        dAdp(6,1) = p.T_f_2_over_h_f;
    elseif strcmp(param_name,'D_d_min')
        dbdp(1) = -1;
    elseif strcmp(param_name,'h')
        dbdp(5:6) = 1;
    elseif strcmp(param_name,'h_d_over_D_s')
        dAdp(7,1) = -1;
    end

end