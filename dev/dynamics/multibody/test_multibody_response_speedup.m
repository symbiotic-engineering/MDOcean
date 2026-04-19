function test_multibody_response_speedup()
% TEST_MULTIBODY_RESPONSE_SPEEDUP  Regression test and benchmark for
% multibody_response with and without the get_multibody_helper_terms fast
% path.
%
% Generates deterministic 14x15 positive-real inputs, calls
% multibody_response twice (standard slow path and pre-computed helper
% path), verifies the outputs match to within floating-point tolerance, and
% reports timeit runtimes for both call patterns.

% Locate production code relative to this file
here   = fileparts(mfilename('fullpath'));
dyndir = fullfile(here, '..', '..', '..', 'mdocean', 'simulation', 'modules', 'dynamics');
addpath(dyndir);

% --- Deterministic 14x15 positive-real inputs ---
rng(0);
sz = [14, 15];
B_c       = rand(sz);
B_f       = rand(sz);
B_s       = rand(sz);
K_f       = rand(sz);
K_s       = rand(sz);
m_c       = rand(sz);
m_f       = rand(sz);
m_s       = rand(sz);
w         = rand(sz);
K_p       = rand(sz);
B_p       = rand(sz);
F_f_mag   = rand(sz);
F_f_phase = rand(sz);
F_s_mag   = rand(sz);
F_s_phase = rand(sz);

% =========================================================================
% Reference: standard call (all 10 outputs, no helper terms)
% =========================================================================
[mag_U_ref, real_P_ref, mag_X_u_ref, mag_X_f_ref, mag_X_s_ref, ...
 phase_X_f_ref, phase_X_s_ref, phase_U_ref, imag_P_ref, phase_X_u_ref] = ...
    multibody_response(B_c, B_f, B_s, K_f, K_s, m_c, m_f, m_s, w, ...
                       K_p, B_p, F_f_mag, F_f_phase, F_s_mag, F_s_phase);

% 3-output reference (tests the nargout guard)
[mag_U_3, real_P_3, mag_X_u_3] = ...
    multibody_response(B_c, B_f, B_s, K_f, K_s, m_c, m_f, m_s, w, ...
                       K_p, B_p, F_f_mag, F_f_phase, F_s_mag, F_s_phase);

% =========================================================================
% Fast path: pre-compute helper terms then call with individual args
% =========================================================================
[h_t18,h_t9,h_t108,h_t121,h_t_110_118,h_t144,h_t109,h_t119,h_D_sys, ...
 h_t145,h_t103,h_t126,h_t130,h_t129,h_t131,h_t102,h_t124] = ...
    get_multibody_helper_terms(B_c, B_f, B_s, K_f, K_s, m_c, m_f, m_s, w, ...
                                F_f_mag, F_f_phase, F_s_mag, F_s_phase);

[mag_U_h, real_P_h, mag_X_u_h, mag_X_f_h, mag_X_s_h, ...
 phase_X_f_h, phase_X_s_h, phase_U_h, imag_P_h, phase_X_u_h] = ...
    multibody_response(B_c, B_f, B_s, K_f, K_s, m_c, m_f, m_s, w, ...
                       K_p, B_p, F_f_mag, F_f_phase, F_s_mag, F_s_phase, [], ...
                       h_t18,h_t9,h_t108,h_t121,h_t_110_118,h_t144,h_t109,h_t119, ...
                       h_D_sys,h_t145,h_t103,h_t126,h_t130,h_t129,h_t131,h_t102,h_t124);

% 3-output fast path
[mag_U_h3, real_P_h3, mag_X_u_h3] = ...
    multibody_response(B_c, B_f, B_s, K_f, K_s, m_c, m_f, m_s, w, ...
                       K_p, B_p, F_f_mag, F_f_phase, F_s_mag, F_s_phase, [], ...
                       h_t18,h_t9,h_t108,h_t121,h_t_110_118,h_t144,h_t109,h_t119, ...
                       h_D_sys,h_t145,h_t103,h_t126,h_t130,h_t129,h_t131,h_t102,h_t124);

% =========================================================================
% Correctness checks
% =========================================================================
tol = 1e-10;

check('mag_U    (10 out): slow vs fast',  mag_U_ref,       mag_U_h,       tol);
check('real_P   (10 out): slow vs fast',  real_P_ref,      real_P_h,      tol);
check('mag_X_u  (10 out): slow vs fast',  mag_X_u_ref,     mag_X_u_h,     tol);
check('mag_X_f  (10 out): slow vs fast',  mag_X_f_ref,     mag_X_f_h,     tol);
check('mag_X_s  (10 out): slow vs fast',  mag_X_s_ref,     mag_X_s_h,     tol);
check('phase_X_f(10 out): slow vs fast',  phase_X_f_ref,   phase_X_f_h,   tol);
check('phase_X_s(10 out): slow vs fast',  phase_X_s_ref,   phase_X_s_h,   tol);
check('phase_U  (10 out): slow vs fast',  phase_U_ref,     phase_U_h,     tol);
check('imag_P   (10 out): slow vs fast',  imag_P_ref,      imag_P_h,      tol);
check('phase_X_u(10 out): slow vs fast',  phase_X_u_ref,   phase_X_u_h,   tol);

check('mag_U   (3 out): slow 3 vs 10',   mag_U_ref,       mag_U_3,       tol);
check('real_P  (3 out): slow 3 vs 10',   real_P_ref,      real_P_3,      tol);
check('mag_X_u (3 out): slow 3 vs 10',   mag_X_u_ref,     mag_X_u_3,     tol);

check('mag_U   (3 out): fast 3 vs 10',   mag_U_h,         mag_U_h3,      tol);
check('real_P  (3 out): fast 3 vs 10',   real_P_h,        real_P_h3,     tol);
check('mag_X_u (3 out): fast 3 vs 10',   mag_X_u_h,       mag_X_u_h3,    tol);

fprintf('All checks passed.\n\n');

% =========================================================================
% Benchmarks
% =========================================================================
slow_fn = @() multibody_response(B_c, B_f, B_s, K_f, K_s, m_c, m_f, m_s, w, ...
                                  K_p, B_p, F_f_mag, F_f_phase, F_s_mag, F_s_phase);

fast_fn = @() multibody_response(B_c, B_f, B_s, K_f, K_s, m_c, m_f, m_s, w, ...
                                  K_p, B_p, F_f_mag, F_f_phase, F_s_mag, F_s_phase, [], ...
                                  h_t18,h_t9,h_t108,h_t121,h_t_110_118,h_t144,h_t109,h_t119, ...
                                  h_D_sys,h_t145,h_t103,h_t126,h_t130,h_t129,h_t131,h_t102,h_t124);

t_slow = timeit(slow_fn);
t_fast = timeit(fast_fn);

fprintf('multibody_response (slow path, no helper): %.3f ms\n', t_slow*1e3);
fprintf('multibody_response (fast path, with h):    %.3f ms\n', t_fast*1e3);
fprintf('Speedup: %.1fx\n', t_slow/t_fast);

end

% -------------------------------------------------------------------------
function check(label, A, B, tol)
    err = max(abs(A(:) - B(:))) ./ (max(abs(A(:))) + eps);
    if err > tol
        error('FAIL: %s  (rel err = %.2e, tol = %.2e)', label, err, tol);
    end
    fprintf('  ok  %s\n', label);
end
