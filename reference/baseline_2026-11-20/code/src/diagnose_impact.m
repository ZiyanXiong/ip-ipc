function D = diagnose_impact(S, P, Rip, Ripc, f, label)
%DIAGNOSE_IMPACT  Same-state replay of frame f and the step-limit analysis at the first limiting IP iterate.
% Both solvers restart from IP's saved frame-start state (positions, velocities, inertial prediction) and IP keeps
% its actual multiplier history. Diagnostic solves (not simulation work) are counted in D.diagnostic_solves.
if nargin < 6, label = 'selected'; end
st = Rip.starts{f};
D = struct('label', label, 'frame', f, 'diagnostic_solves', 0);
Fw = S.F_w; dhat = S.d_hat; tol = S.tol_disp; gy = S.ground_y;

% ---- 1. replays from the identical state
[~, Fi, iti] = ip_solve_frame(S, P, st.x, st.xt, f, st.hist, true);
[~, Fc, itc] = ipc_solve_frame(S, P, st.x, st.xt, f, true);
ref = Rip.iters([Rip.iters.frame] == f);
same_len = numel(iti) == numel(ref);
D.ip_replay = struct('main', Fi.main, 'verification', Fi.verification, 'status', Fi.status, ...
    'trajectory_main', Rip.frames(f).main, 'trajectory_verification', Rip.frames(f).verification, ...
    'same_iteration_count', same_len, ...
    'max_abs_diff_alpha', ifelse(same_len, @() max(abs(nan0([iti.alpha]) - nan0([ref.alpha]))), NaN), ...
    'max_abs_diff_rd', ifelse(same_len, @() max(abs(nan0([iti.rd_inf]) - nan0([ref.rd_inf]))), NaN));
D.ipc_replay = struct('main', Fc.main, 'status', Fc.status, 'trajectory_main_on_own_trajectory', Ripc.frames(f).main);
% identical entry state: both replays' first iterates and the stored prediction
D.same_state = struct('max_abs_diff_first_iterate_x', max(abs(iti(1).cap.x - itc(1).cap.x)), ...
    'max_abs_diff_xt_vs_prediction', max(abs(iti(1).cap.xt - (st.x + S.dt * st.v + S.dt2 * S.gvec))), ...
    'hist_rows', numel(field_or(st.hist, 'rv', [])), 'hist_frame', field_or(st.hist, 'frame', NaN), ...
    'rows_carried_at_start', Fi.hist_reused, 'rows_central_at_start', Fi.hist_fresh);
% supplementary: the same frame with every multiplier centrally started (history ignored), identical state
[~, Fz] = ip_solve_frame(S, P, st.x, st.xt, f, [], false);
D.ip_central_start_replay = struct('main', Fz.main, 'verification', Fz.verification, 'status', Fz.status);
D.replay_solves = Fi.main + Fi.verification + Fc.main + Fz.main + Fz.verification;
D.ip_iters = strip_cap(iti); D.ipc_iters = strip_cap(itc);

% ---- 2. the first limiting IP iterate
mains = iti(strcmp({iti.kind}, 'main'));
q = find([mains.alpha_eff] < 1 & [mains.alpha_eff] > 0, 1);
if isempty(q)
    D.limiting = 'none: every accepted step of this frame is the full step'; return;
end
it = mains(q); c = it.cap;
D.iterate = struct('k', it.k, 'solve', it.solve, 'tau_hat', it.tau_hat, 'rows', numel(c.rv), ...
    'accepted_alpha', it.alpha, 'accepted_alpha_eff', it.alpha_eff, 'limiter_recorded', it.limiter, ...
    'levels', c.levels);
x = c.x; xt = c.xt; rv = c.rv; s = c.s; lam = c.lam; tau = c.tau; g = c.g; mu = tau;
y = x(S.ydof) - gy;
[E, gE, H0] = mech_eval(S, x, xt, 3);
JT = @(r) accum_y(S, rv, r);
w = lam ./ s;
K = H0 + sparse(2*rv, 2*rv, w, S.ndof, S.ndof);
b = gE - JT(mu ./ s) + JT(w .* (g - s));
dx = shared_solve(K, -b); D.diagnostic_solves = D.diagnostic_solves + 1;
D.recomputed_direction_max_abs_diff = max(abs(dx - c.dx));
ds = dx(2*rv) + (g - s);
dlam = mu ./ s - lam - w .* ds;

% ---- 3. assembly checks (s = g at this iterate?) and the unreduced KKT solve
A = struct();
A.max_abs_g_minus_s_dhat = max(abs(g - s)) / dhat;
dK = K - H0;
A.stiffness_rel_err = max(abs(full(diag(dK(2*rv, 2*rv))) - lam ./ s) ./ (lam ./ s));
dK(sub2ind(size(dK), 2*rv, 2*rv)) = 0;
A.offdiag_contact_terms = full(max(abs(dK(:))));
db = b - gE;
A.rhs_effective_force_rel_err = max(abs(-db(2*rv) - mu ./ s) ./ (mu ./ s));   % valid when s = g
A.rhs_uses_actual_lambda_would_differ_by_Fw = max(abs(mu ./ s - lam)) / Fw;
rd = gE - JT(lam);
A.rd_inf_actual_lambda_Fw = max(abs(rd)) / Fw;
A.rd_inf_effective_force_Fw = max(abs(gE - JT(mu ./ s))) / Fw;
% full Newton system on (dx, ds, dlambda): [H0 0 -J'; J -I 0; 0 Lam S] = -[rd; g - s; s lam - mu]
m = numel(rv); n = S.ndof;
Jm = sparse(1:m, 2*rv, 1, m, n);
KKT = [H0, sparse(n, m), -Jm'; Jm, -speye(m), sparse(m, m); sparse(m, n), spdiags(lam, 0, m, m), spdiags(s, 0, m, m)];
rhs = -[rd; g - s; s .* lam - mu];
z = KKT \ rhs; D.diagnostic_solves = D.diagnostic_solves + 1;
A.kkt_dx_rel_err = norm(z(1:n) - dx) / norm(dx);
A.kkt_ds_rel_err = norm(z(n+1:n+m) - ds) / max(norm(ds), eps);
A.kkt_dlambda_rel_err = norm(z(n+m+1:end) - dlam) / max(norm(dlam), eps);
D.assembly = A;

% ---- 4. the limiting vertex and every restriction, from the raw direction
eta = P.ip.ftb_eta;
R = step_limits(S, P, x, xt, rv, s, lam, tau, g, y, dx, ds, dlam, c.nu, E, gE, it.split);
j = R.slack_row;
lv = rv(j);
[~, ~, ~, Cipc] = ipc_contact_eval(S, x, 0);
kn_ipc = NaN; fn_ipc = NaN;
ia = find(Cipc.act == lv); if ~isempty(ia), kn_ipc = Cipc.kn(ia); fn_ipc = Cipc.fn(ia); end
V = struct();
V.vertex = lv; V.row = j; V.is_bottom_vertex = ismember(lv, S.bottom);
V.source = src_name(c.src(j));
V.g_dhat = g(j) / dhat; V.s_dhat = s(j) / dhat;
V.lambda_Fw = lam(j) / Fw; V.mu_FwL = mu / (Fw * S.L);
V.s_lambda_over_mu = s(j) * lam(j) / mu;
V.lambda_over_s_dhat_over_Fw = lam(j) / s(j) * dhat / Fw;           % condensed normal stiffness
V.lambda_over_s_over_mv = lam(j) / s(j) / S.m(lv);                   % relative to the vertex's inertia block
V.central_stiffness_mu_over_s2_dhat_over_Fw = mu / s(j)^2 * dhat / Fw;
V.stiffness_ratio_lambda_s_to_central = (lam(j) / s(j)) / (mu / s(j)^2);
V.effective_force_mu_over_s_Fw = mu / s(j) / Fw;
V.ipc_force_same_gap_Fw = -2 * S.kappa_dt2 * barrier_dB(g(j), dhat) * g(j) / Fw;
V.ipc_stiffness_same_gap_dhat_over_Fw = ipc_kn(S, g(j)) * dhat / Fw;
V.Jdx_dhat = dx(2*lv) / dhat;
V.predicted_motion_over_gap = dx(2*lv) / g(j);   % negative: closing
V.ds_dhat = ds(j) / dhat; V.dlambda_Fw = dlam(j) / Fw;
V.raw_slack_crossing = s(j) / (-ds(j));
V.slack_bound = eta * s(j) / (-ds(j));
% the non-contact push on that vertex (along n) and its non-contact stiffness
V.noncontact_force_Fw = -gE(2*lv) / Fw;
V.noncontact_stiffness_over_mv = full(H0(2*lv, 2*lv)) / S.m(lv);
D.vertex = V; D.restrictions = R;

% ---- 5. the diagnostic direction: identical positions, rows, slacks and mu; lambda_i := mu / s_i
lam_c = mu ./ s; w_c = lam_c ./ s;
K_c = H0 + sparse(2*rv, 2*rv, w_c, S.ndof, S.ndof);
b_c = gE - JT(mu ./ s) + JT(w_c .* (g - s));
dx_c = shared_solve(K_c, -b_c); D.diagnostic_solves = D.diagnostic_solves + 1;
ds_c = dx_c(2*rv) + (g - s);
dlam_c = mu ./ s - lam_c - w_c .* ds_c;
nu_c = P.ip.nu_rho * max(abs(lam_c + dlam_c));
Rc = step_limits(S, P, x, xt, rv, s, lam_c, tau, g, y, dx_c, ds_c, dlam_c, nu_c, E, gE, it.split);
C = struct();
C.same_vertex_predicted_motion_over_gap = dx_c(2*lv) / g(j);
C.same_vertex_Jdx_dhat = dx_c(2*lv) / dhat;
C.same_vertex_stiffness_lambda_over_s_dhat_over_Fw = w_c(j) * dhat / Fw;
C.limits = Rc;
C.limiting_vertex = NaN; if Rc.slack_row > 0, C.limiting_vertex = rv(Rc.slack_row); end
C.max_closure_over_gap = max(-dx_c(2*rv) ./ g);
C.direction_rel_change = norm(dx_c - dx) / norm(dx);
C.rd_inf_with_centered_lambda_Fw = max(abs(gE - JT(lam_c))) / Fw;
D.centered = C;

% ---- 6. IPC at the same geometry (its own active set, force and stiffness at x)
[Em, gm, Hm] = mech_eval(S, x, xt, 3);
[Eb, gb, Hb, Cb] = ipc_contact_eval(S, x, 3);
dx_i = shared_solve(Hm + Hb, -(gm + gb)); D.diagnostic_solves = D.diagnostic_solves + 1;
I = struct();
I.active = numel(Cb.act); I.same_pairs_as_ip_rows = isequal(sort(Cb.act), sort(rv));
I.rd_inf_Fw = max(abs(gm + gb)) / Fw; I.dx_tol = max(abs(dx_i)) / tol;
I.same_vertex_predicted_motion_over_gap = dx_i(2*lv) / g(j);
I.same_vertex_Jdx_dhat = dx_i(2*lv) / dhat;
I.same_vertex_stiffness_dhat_over_Fw = kn_ipc * dhat / Fw; I.same_vertex_force_Fw = fn_ipc / Fw;
dyi = dx_i(S.ydof); idx = find(dyi < 0);
[tm, qq] = min(y(idx) ./ (-dyi(idx)) * (1 - P.ipc.ccd_eta));
if isempty(tm) || tm >= 1.1, I.ccd = 1; I.ccd_vertex = NaN; else, I.ccd = min(1, tm); I.ccd_vertex = idx(qq); end
fe = elem_check(S, x, dx_i, I.ccd);
I.elem = Inf; if fe <= I.ccd, I.elem = P.valid_eta * fe; end
a = min(I.ccd, I.elem);
E0 = Em + Eb; I.trials = [];
for t = 1:P.ls_max
    xtr = x + a * dx_i;
    [~, np] = elem_check(S, xtr);
    if np > 0, I.trials = [I.trials; a, NaN, 0]; a = a / 2; continue; end
    Et = mech_eval(S, xtr, xt, 1) + ipc_contact_eval(S, xtr, 1);
    I.trials = [I.trials; a, Et - E0, Et <= E0];
    if Et <= E0, break; end
    a = a / 2;
end
I.accepted_alpha = a;
I.max_closure_over_gap_ip_rows = max(-dx_i(2*rv) ./ g);
D.ipc_same_geometry = I;

% ---- 7. per-row table (carried, centered, IPC)
T = table(rv, g / dhat, s / dhat, lam / Fw, s .* lam / mu, w * dhat / Fw, (mu ./ s.^2) * dhat / Fw, ...
    dx(2*rv) ./ g, dx_c(2*rv) ./ g, dx_i(2*rv) ./ g, eta * s ./ max(-ds, 0), eta * s ./ max(-ds_c, 0), ...
    arrayfun(@(v) ipc_kn(S, y(v)), rv) * dhat / Fw, arrayfun(@(q2) string(src_name(c.src(q2))), (1:numel(rv))'), ...
    'VariableNames', {'vertex', 'g_dhat', 's_dhat', 'lambda_Fw', 's_lambda_over_mu', 'lambda_over_s_dhat_Fw', ...
    'mu_over_s2_dhat_Fw', 'motion_over_gap_carried', 'motion_over_gap_centered', 'motion_over_gap_ipc', ...
    'slack_bound_carried', 'slack_bound_centered', 'ipc_kn_dhat_Fw', 'source'});
D.rows_table = T;
end

% ======================================================================================== helpers
function R = step_limits(S, P, x, xt, rv, s, lam, tau, g, y, dx, ds, dlam, nu, E, gE, split_recorded)
eta = P.ip.ftb_eta; dhat = S.d_hat;
R = struct();
[R.slack, R.slack_row] = ftbm(s, ds, eta);
[R.multiplier, R.multiplier_row] = ftbm(lam, dlam, eta);
R.slack = min(1, R.slack); R.multiplier = min(1, R.multiplier);
R.slack_vertex = NaN; if R.slack_row > 0, R.slack_vertex = rv(R.slack_row); end
R.multiplier_vertex = NaN; if R.multiplier_row > 0, R.multiplier_vertex = rv(R.multiplier_row); end
dy = dx(S.ydof); idx = find(dy < 0);
[tg, qg] = min(y(idx) ./ (-dy(idx)));
if isempty(tg), R.ground_exact = Inf; R.ground_vertex = NaN; else, R.ground_exact = tg; R.ground_vertex = idx(qg); end
R.ground_vertex_has_row = ismember(R.ground_vertex, rv);
fe = elem_check(S, x, dx, 1);
R.element = Inf; if fe <= 1, R.element = P.valid_eta * fe; end
% the proposal as the solver forms it
levels0 = tau > S.tau_final * (1 + 1e-12);   % first level (tau_hat = 1)
a = R.slack; if R.multiplier < a, a = R.multiplier; end
split_now = P.ip.split_tau1 && levels0 && R.multiplier < R.slack;
if split_now, a = R.slack; end
fe2 = elem_check(S, x, dx, a);
if fe2 <= a, a = P.valid_eta * fe2; end
R.proposal = a; R.split = split_now;
% merit backtracking from the proposal
rp = g - s; rp1 = sum(abs(rp));
slope = gE' * dx - tau * sum(ds ./ s) - nu * rp1;
phi0 = E - tau * sum(log(s / dhat)) + nu * rp1;
R.trials = [];
acc = NaN;
for t = 1:P.ls_max
    xtr = x + a * dx;
    [~, np] = elem_check(S, xtr);
    if np > 0, R.trials = [R.trials; a, NaN, NaN, 0]; a = a / 2; continue; end
    Et = mech_eval(S, xtr, xt, 1);
    gtr = xtr(2*rv) - S.ground_y;
    str = max(s + a * ds, gtr);
    phi = Et - tau * sum(log(str / dhat)) + nu * sum(abs(gtr - str));
    ok = phi <= phi0 + P.ip.armijo_c * a * slope;
    R.trials = [R.trials; a, phi - phi0, P.ip.armijo_c * a * slope, ok];
    if ok, acc = a; break; end
    a = a / 2;
end
R.merit_accepted = acc;
R.merit_backtracks = size(R.trials, 1) - 1;
% the early ground scan at the accepted point
if isfinite(acc)
    yacc = y + acc * dy;
    R.scan_all_positive = all(yacc > 0);
    R.min_gap_after_dhat = min(yacc) / dhat;
else
    R.scan_all_positive = NaN; R.min_gap_after_dhat = NaN;
end
% which restriction set the accepted displacement. IP applies slack, multiplier (update B only; not under the
% split rule) and element bounds a priori, then merit backtracking; the all-vertex ground safeguard acts after the
% accept (fallback halving from the last checked position) and is reported separately.
cands = [R.slack, R.multiplier * ~split_now + split_now, R.element];
names = {'slack positivity', 'multiplier positivity', 'element validity'};
[~, qm] = min(cands);
if isfinite(acc) && ~R.scan_all_positive
    R.binding = 'ground safeguard (fallback after the accept)';
elseif isfinite(acc) && acc < R.proposal
    R.binding = 'merit backtracking';
elseif R.proposal >= 1
    R.binding = 'none (full step)';
else
    R.binding = names{qm};
end
R.split_recorded = split_recorded;
end

function [a, j] = ftbm(v, dv, eta)
a = Inf; j = 0;
idx = find(dv < 0);
if isempty(idx), return; end
[a, q] = min(-eta * v(idx) ./ dv(idx)); j = idx(q);
end

function v = accum_y(S, rv, r)
v = zeros(S.ndof, 1); v(2*rv) = r;
end

function d = barrier_dB(gap, dhat)
[~, d] = ipc_barrier(gap^2, dhat^2);
end

function k = ipc_kn(S, gap)
if gap >= S.d_hat || gap <= 0, k = 0; return; end
[~, dB, ddB] = ipc_barrier(gap^2, S.d_hat^2);
k = S.kappa_dt2 * (4 * gap^2 * ddB + 2 * dB);
end

function n = src_name(code)
names = {'central (frame start)', 'carried (history)', 'early insertion (central)', 'outer refresh (central)'};
n = names{code + 1};
end

function v = nan0(v)
v(isnan(v)) = 0;
end

function v = ifelse(c, f, d)
if c, v = f(); else, v = d; end
end

function its = strip_cap(its)
for q = 1:numel(its), its(q).cap = []; end
end

function v = field_or(s, f, d)
if isempty(s) || ~isfield(s, f), v = d; else, v = s.(f); end
end
