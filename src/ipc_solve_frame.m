function [x, F, its] = ipc_solve_frame(S, P, xn, xt, frame, capture)
%IPC_SOLVE_FRAME  One implicit-Euler frame of native IPC (frictionless ground) with the matched stop.
% Newton on E_mech + kappa dt^2 sum B(D) with the active ground pairs evaluated at every iterate and trial.
% Step: full-step half-plane CCD (0.9 of each closing vertex's gap; CFL off; the simplex feasible-step cap has no
% pairs here), the shared element-validity bound, then backtracking on the incremental energy (E <= E0, halving,
% line_search/max_iter trials; on exhaustion the last tried trial is kept, as natively).
% Matched stop (round 43 'ipc'): the native displacement test |dx|_inf <= tol_disp AND, at that same state,
% ||r_d||_inf <= residual_target F_w with IPC's own law, full-correction CCD >= ccd_tol, valid elements and
% positive gaps. The state is returned without the terminal step, so its correction is the one checked.
if nargin < 6, capture = false; end
x = xn;
F = struct('solver', 'ipc', 'frame', frame, 'status', 'running', 'reason', '', 'main', 0, 'verification', 0, ...
    'discarded', 0, 'fallback', 0, 'ls_limit', 0, 'cert_rejects', 0, 'matched_rejects', 0, 'matched_reject_rd', 0, ...
    'matched_reject_ccd', 0, 'matched_reject_geo', 0, 'trials', 0);
its = [];
tol = S.tol_disp; Fw = S.F_w;
for k = 0:P.cap_main-1
    [Em, gm, Hm] = mech_eval(S, x, xt, 3);
    [Eb, gb, Hb, C] = ipc_contact_eval(S, x, 3);
    grad = gm + gb; K = Hm + Hb;
    dx = shared_solve(K, -grad);
    F.main = F.main + 1;
    dx_inf = max(abs(dx));
    rd_inf = max(abs(grad)) / Fw;
    y = x(S.ydof) - S.ground_y; dy = dx(S.ydof);
    it = it_new(frame, k, F.main);
    it.rows = numel(C.act); it.rd_inf = rd_inf; it.dx_tol = dx_inf / tol;
    it.min_gap_dhat = min(y) / S.d_hat; it.E = Em + Eb;
    if capture
        it.cap = struct('x', x, 'dx', dx, 'act', C.act, 'fn', C.fn, 'kn', C.kn);
    end
    if dx_inf <= tol
        [full_toi, ~] = ccd_toi(y, dy, P.ipc.ccd_eta);
        [~, np] = elem_check(S, x);
        geo = np == 0 && all(y > 0);
        first_full = elem_check(S, x, dx, 1);
        it.native_stop = true; it.full_ccd = full_toi; it.full_valid = isinf(first_full);
        pass = rd_inf <= P.residual_target && full_toi >= P.ipc.ccd_tol && geo;
        if pass
            it.event = 'matched stop: returned';
            it.alpha = 0; it.limiter = 'stop';
            its = [its, it]; %#ok<AGROW>
            F.status = 'converged'; F.reason = 'matched stop';
            break;
        end
        F.matched_rejects = F.matched_rejects + 1;
        F.matched_reject_rd = F.matched_reject_rd + (rd_inf > P.residual_target);
        F.matched_reject_ccd = F.matched_reject_ccd + (full_toi < P.ipc.ccd_tol);
        F.matched_reject_geo = F.matched_reject_geo + ~geo;
        it.event = 'native stop rejected by the matched check';
    end
    [a_ccd, vccd] = ccd_toi(y, dy, P.ipc.ccd_eta);
    a_ccd = min(1, a_ccd);
    alpha = a_ccd;
    first = elem_check(S, x, dx, alpha);
    a_elem = Inf;
    if first <= alpha
        a_elem = P.valid_eta * first; alpha = a_elem;
    end
    a_prop = alpha;
    E0 = Em + Eb;
    accepted = false; last_ok = NaN; cert = 0; tr = 0; Etr = NaN;
    for tr = 1:P.ls_max
        xtr = x + alpha * dx;
        [~, np] = elem_check(S, xtr);
        if np > 0
            cert = cert + 1; alpha = alpha / 2; continue;
        end
        Etr = mech_eval(S, xtr, xt, 1) + ipc_contact_eval(S, xtr, 1);
        last_ok = alpha;
        if Etr <= E0
            accepted = true; break;
        end
        alpha = alpha / 2;
    end
    F.trials = F.trials + tr; F.cert_rejects = F.cert_rejects + cert;
    if ~accepted
        F.ls_limit = F.ls_limit + 1;
        if ~isnan(last_ok)
            alpha = last_ok; x = x + alpha * dx;  % native: the last tried trial is kept (a warning)
        else
            alpha = 0;
            F.discarded = F.discarded + 1;
        end
        it.event = strjoin({it.event, 'line search limit'}, '; ');
    else
        x = xtr;
    end
    it.alpha_ccd = a_ccd; it.ccd_vertex = vccd; it.alpha_elem = a_elem; it.alpha_proposal = a_prop;
    it.alpha = alpha; it.alpha_eff = alpha; it.trials = tr; it.E_trial = Etr;
    if alpha >= 1
        it.limiter = 'none';
    elseif alpha < a_prop
        it.limiter = 'energy backtracking';
    elseif a_elem <= a_ccd
        it.limiter = 'element validity';
    else
        it.limiter = 'ccd';
    end
    its = [its, it]; %#ok<AGROW>
end
if strcmp(F.status, 'running')
    F.status = 'budget'; F.reason = 'main-system budget reached';
end
end

function it = it_new(frame, k, solve)
it = struct('solver', 'ipc', 'frame', frame, 'k', k, 'solve', solve, 'kind', 'main', 'tau_hat', NaN, 'rows', 0, ...
    'rd_inf', NaN, 'rp_inf_dhat', NaN, 'cent_hat', NaN, 'dx_tol', NaN, 'alpha_slack', NaN, 'alpha_mult', NaN, ...
    'alpha_ccd', NaN, 'alpha_elem', NaN, 'alpha_proposal', NaN, 'alpha', NaN, 'alpha_lambda', NaN, 'fallback_t', NaN, ...
    'alpha_eff', NaN, 'trials', 0, 'limiter', '', 'lim_vertex', NaN, 'ccd_vertex', NaN, 'split', false, ...
    'native_stop', false, 'full_ccd', NaN, 'full_valid', NaN, 'min_gap_dhat', NaN, 'E', NaN, 'E_trial', NaN, ...
    'event', '', 'inserted', 0, 'cap', []);
end

function [toi, v] = ccd_toi(y, dy, eta)
% easy_vertex_half_plane_trajectory_filter.cu: toi = g0 / (-n.dx) (1 - eta) for every vertex moving toward the
% plane; 1.1 if none
idx = find(dy < 0);
toi = 1.1; v = NaN;
if isempty(idx), return; end
t = y(idx) ./ (-dy(idx)) * (1 - eta);
[tm, j] = min(t);
if tm < toi
    toi = tm; v = idx(j);
end
end
