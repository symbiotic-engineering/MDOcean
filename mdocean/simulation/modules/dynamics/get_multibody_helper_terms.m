function [t18,t9,t108,t121,t_110_118,t144,t109,t119,h_D_sys,t145,t103,t126,t130,t129,t131,t102,t124] = ...
    get_multibody_helper_terms(B_c,B_f,B_s,K_f,K_s,m_c,m_f,m_s,w, ...
                               F_f_mag,F_f_phase,F_s_mag,F_s_phase,D_sys)
% GET_MULTIBODY_HELPER_TERMS  Pre-compute all controller-independent
% intermediate variables for MULTIBODY_RESPONSE.
%
% The 17 return values may be passed as args 17–33 of MULTIBODY_RESPONSE so
% that these quantities are not recomputed on every controller evaluation in
% the brute-force optimisation loop.  K_p and B_p are the only inputs of
% MULTIBODY_RESPONSE that are absent here; all calculations below are
% therefore independent of the controller.
%
% D_sys is optional: when omitted it is computed internally and returned as
% h_D_sys (the 9th output).
%
% Usage:
%   [t18,t9,t108,t121,t_110_118,t144,t109,t119,h_D_sys, ...
%    t145,t103,t126,t130,t129,t131,t102,t124] = ...
%       get_multibody_helper_terms(B_c,B_f,B_s,K_f,K_s,m_c,m_f,m_s,w, ...
%                                  F_f_mag,F_f_phase,F_s_mag,F_s_phase);
%   [mag_U,real_P,mag_X_u] = multibody_response(B_c,B_f,B_s,K_f,K_s, ...
%                               m_c,m_f,m_s,w,K_p,B_p, ...
%                               F_f_mag,F_f_phase,F_s_mag,F_s_phase,[], ...
%                               t18,t9,t108,t121,t_110_118,t144,t109,t119, ...
%                               h_D_sys,t145,t103,t126,t130,t129,t131,t102,t124);

% ---- variables that appear before K_p is first used in multibody_response ----
s = w*1i;
t3 = B_c.*w;
t4 = B_f.*w;
t5 = B_s.*w;
K_s_B_f = K_s .* B_f;
K_f_B_s = K_f .* B_s;
t6  = K_f.*w;
t8  = m_c.^2;
t9  = w.^2;
t10 = t9 .* w;
t12 = K_f.*K_s;
t16 = -K_f;
t17 = -K_s;
t18 = 1.0./w;
t30 = complex(cos(F_f_phase), sin(F_f_phase));
t31 = complex(cos(F_s_phase), sin(F_s_phase));
t27 = K_f.*1i;
t28 = K_s.*1i;
t34 = m_c.*s;
t11 = t9.^2;
t13 = t3.*2.0;
t14 = K_s_B_f .* w;
t15 = K_f_B_s .* w;
t19 = -t6;
t21 = m_c.*t10;
t22 = m_f.*t9;
t23 = m_f.*t10;
t24 = m_s.*t9;
t29 = -t12;
t33 = t5.*1i;
t37 = m_c.*t9.*2.0;
t40 = -t27;
t41 = -t28;
t42 = t4.*t5;
t45 = m_s.*t6.*w;
t47 = t12.*1i;
t56 = t3.^2;
t62 = t3.*s;
t63 = t4.*s;
t65 = 1i*t37;
t64 = t65/2;
t68 = B_c+t34;

% ---- remaining controller-independent terms ----
t32 = m_f.*m_s.*t11;
t43 = t4.*t24;
t44 = t5.*t22;
t46 = K_s.*t22;
t48 = t8.*t11;
t49 = 1.0./t30;
t50 = 1.0./t31;
t51 = K_s_B_f .* s;
t52 = K_f_B_s .* s;
t53 = -t21;
t54 = -t23;
t55 = -t24;
t57 = m_c.*t9.*t13;
t58 = -t47;
t66 = t22.*1i;
t67 = t24.*1i;
t74 = -t64;
t75 = -t65;
t78 = t4.*t33;
t79 = -t56;
t80 = t3.*t65;
t83 = t45.*1i;
t84 = t22.*t28;
t88 = t56.*1i;
t93 = t4+t5+t13;
t94 = -t88;
t95 = t3+t64;
t100 = t16+t17+t22+t24+t37;
t105 = t19+t21+t23+t62+t63;
t59 = -t32;
t70 = -t43;
t71 = -t44;
t72 = t32.*1i;
t76 = -t66;
t77 = -t67;
t81 = t43.*1i;
t85 = t48.*1i;
t91 = t44.*-1i;
t86 = -t72;
t90 = -t81;

if nargin < 14 || isempty(D_sys)
    D_sys = t12+t32-t42-t48+t51+t52+t56+t80+t90+t91+t17.*t22+m_s.*t19.*w;
end

t97  = t93.^2;
t98  = K_s+t33+t55;
t99  = t4+t40+t66;
t101 = t100.^2;
t103 = t5+t41+t67+t95;
t106 = t6+t53+t54+t62+t63;
t109 = t40+t41+t65+t66+t67+t93;
t102 = t95+t99;
t104 = t3+t5+t28+t74+t77;
t107 = t97+t101;
t110 = t27+t28+t75+t76+t77+t93;

t111 = 1.0./t109;
t108 = sqrt(t107);
t112 = conj(t111);
t_Ff30    = F_f_mag.*t30;
t113      = t_Ff30.*t103.*t111;
t115      = t29+t42+t45+t46+t48+t51+t52+t59+t79+t80+t90+t91;
t117      = abs(t115);
t119      = -1.0./D_sys;
t122      = t14+t15+t57+t58+t70+t71+t78+t83+t84+t85+t86+t94;
t114      = F_f_mag.*t49.*t104.*t112;
t118      = 1.0./t115;
t121      = 1.0./t117;
t124      = 1.0./t122;
t_110_118 = t110.*t118;
t_Fs31_124   = F_s_mag.*t31.*t124;
t_Ff30_124w  = t_Ff30.*t124.*w;
t142      = (F_s_mag.*t18.*t50.*t105.*t112.*t115)./complex(imag(D_sys),real(D_sys));
t_Fs31_t9_t68 = t_Fs31_124.*t9.*t68;
t126      = complex(-imag(t_Fs31_t9_t68), real(t_Fs31_t9_t68));
t127      = t_Ff30_124w.*t98;
t_Ff30w_t95  = t_Ff30_124w.*t95;
t128      = complex(-imag(t_Ff30w_t95), real(t_Ff30w_t95));
t_Fs31_t99w  = t_Fs31_124.*t99.*w;
t131      = complex(-imag(t_Fs31_t99w), real(t_Fs31_t99w));
t143      = t_Fs31_124.*t18.*t106.*t111.*D_sys;
t145      = t114+t142;
t129      = -t128;
t130      = -t127;
t144      = t113+t143;

h_D_sys = D_sys;
